"""Functions used to read input configurations and print trajectories
using the Amber NetCDF format
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

from ipi.utils.units import Constants
import ipi.utils.mathtools as mt

try:
    import scipy
except ImportError:
    scipy = None

__all__ = ["print_nc", "read_nc", 'nc_writer']


def _scipycheck():
    if scipy is None:
        raise RuntimeError(
            """scipy must be installed to use the mode 'nc'. Please run

$~ pip install -U scipy

to install it."""
        )

    # from packaging import version
    #
    # assert version.parse(scipy.__version__) >= version.parse(
    #    "3.22.1"
    # ), "Please get a newer version of scipy"




def NewAmberNetCDF(filename,natoms):
    from scipy.io import netcdf_file
    
    f = netcdf_file(filename,"w")

    f.Conventions = 'AMBER'
    f.ConventionVersion = '1.0'
    f.application = 'AMBER'
    f.program = 'i-pi'
    f.programVersion = 'unknown'
    f.filename = filename
    f.createDimension('frame', None)
    f.createDimension('spatial', 3)
    f.createDimension('cell_spatial', 3)
    f.createDimension('cell_angular', 3)
    f.createDimension('label', 5)
    f.createDimension('atom', natoms)

    time    = f.createVariable('time','f',
                               ('frame',))
    time.units = "picosecond"
    
    spatial = f.createVariable('spatial','c',
                               ('spatial',))
    spatial[:] = ['x','y','z']
    
    coordinates  = f.createVariable('coordinates','f',
                                    ('frame', 'atom', 'spatial'))
    coordinates.units = "angstrom"
    
    cell_spatial = f.createVariable('cell_spatial', 'c',
                                    ('cell_spatial',))
    cell_spatial[:] = ['a','b','c']
    
    cell_angular = f.createVariable('cell_angular', 'c',
                                    ('cell_angular', 'label'))
    cell_angular[:,:] = [['a', 'l', 'p', 'h', 'a'],
                         ['b', 'e', 't', 'a', ' '],
                         ['g', 'a', 'm', 'm', 'a']]

    cell_lengths = f.createVariable('cell_lengths', 'd',
                                    ('frame', 'cell_spatial'))
    cell_lengths.units = "angstrom"
    
    cell_angles  = f.createVariable('cell_angles', 'd',
                                    ('frame', 'cell_angular'))
    cell_angles.units = "degree"
    
    return f


def AppendAmberFrame(netcdfobj,crds,cell,filename,time=None):
    import os
    from scipy.io import netcdf_file
    import numpy as np
    
    if netcdfobj is None:
        if os.path.exists(filename):
            netcdfobj = netcdf_file(filename,"a")
        else:
            netcdfobj = NewAmberNetCDF(filename,crds.shape[0])
    idx = netcdfobj.variables['time'].shape[0]
    if time is None:
        time = idx
    netcdfobj.variables['time'][idx] = time
    netcdfobj.variables['coordinates'][idx,:,:] = crds[:,:].copy()
    netcdfobj.variables['cell_lengths'][idx,:] = np.array([cell[0],cell[1],cell[2]])
    netcdfobj.variables['cell_angles'][idx,:] = np.array([cell[3],cell[4],cell[5]])
    return netcdfobj



from .... engine.outputs import BaseOutput

class NCOutput(BaseOutput):
    
    def __init__(self,filename="out"):
        super(NCOutput, self).__init__(filename=filename)        
        _scipycheck()

    def force_flush(self):
        if self.out is not None:
            self.out.flush()
        
    def flush(self):
        if self.out is not None:
            self.out.flush()
        
    def append(self,crds,cell):
        self.out = AppendAmberFrame(self.out, \
                                    crds,cell,self.filename)
        
    def write(self,data):
        pass

    def fileno(self):
        if hasattr(self.out,'fp'):
            return self.out.fp.fileno()
        else:
            import sys
            return sys.stdout.fileno()

    def close(self):
        self.out.close()

        
        
class NCInput(object):
    
    def __init__(self,filename="out"):
        _scipycheck()
        from scipy.io import netcdf_file
        self.name = filename
        self.out = netcdf_file(filename,"r")
        self.shape = (0,0,0)
        self.cframe = 0
        if 'coordinates' in self.out.variables:
            self.shape = tuple(self.out.variables['coordinates'].shape)

    def close(self):
        self.shape = None
        self.cframe = None
        self.out.close()
        

def print_nc(atoms, cell, filedesc=None, title="", cell_conv=1.0, atoms_conv=1.0):
    """Prints an atomic configuration into an Amber NetCDF trajectory file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable ptraj.TrajectoryWriter object. Defaults to None.
        title: This gives a string to be appended to the comment line.
    """
    filedesc.append( atoms.q.reshape((-1, 3)) * atoms_conv,
                     [ x for x in mt.h2abc_deg(cell.h * cell_conv) ] )
    


def read_nc(filedesc):
    """Reads an Amber NetCDF (.nc) file and returns data in raw format 
    for further units transformation and other post processing. The
    Amber NetCDF format does not store element nor mass information
    (these are stored in a separate "parameter file"), so this routine
    cannot be used to initialize an i-pi simulation.  This routine can
    be used to multiplex several trajectories, however.

    Args:
        filedesc: An open readable pytraj.Trajectory object.

    Returns:
        i-Pi comment line, cell array, data (positions, forces, etc.), atoms names and masses
    """

    import numpy as np

    i = filedesc.cframe
    if i < filedesc.shape[0]:
        comment = "Structure read from Amber NetCDF Format"
        nat = filedesc.shape[1]
        #####################################
        # The topology file is unavailable,
        # so this is not useful as an input
        # format, but it can still be used
        # to multiplex several trajectories
        masses = np.zeros( (nat,) )
        names = ["XX"]*nat
        #####################################
        qatoms = filedesc.out.variables['coordinates'][i,:,:].copy().flatten()
        cell = None
        if 'cell_lengths' in filedesc.out.variables:
            box = np.zeros( (6,) )
            box[0:3] = filedesc.out.variables['cell_lengths'][i,:].copy()
            box[3:6] = [np.deg2rad(x) for x in filedesc.out.variables['cell_angles'][i,:]]
            cell = mt.abc2h(box[0], box[1], box[2], box[3], box[4], box[5])
        filedesc.cframe += 1
        return comment, cell, qatoms, names, masses
    else:
        raise EOFError
    
