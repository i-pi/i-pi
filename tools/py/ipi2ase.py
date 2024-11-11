#!/usr/bin/env python
import argparse
import glob
import os
import re
from warnings import warn
from typing import List, Union, Any, Match
from ase import Atoms
from ase.io import read
from dataclasses import dataclass
from io import TextIOWrapper
import numpy as np
from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.mathtools import abc2h

#---------------------------------------#
# Description of the script's purpose
description = "Convert a i-PI trajectory to a ASE file."

deg2rad     = np.pi / 180.0
abcABC      = re.compile(r"CELL[\(\[\{]abcABC[\)\]\}]: ([-+0-9\.Ee ]*)\s*")
abcABCunits = re.compile(r'\{([^}]+)\}')

#---------------------------------------#
def intype(argument):
    return argument

def itype(x:str):
    return int(x) if x.isdigit() else x


def str2bool(v:Union[bool,str]):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

def is_convertible_to_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False
    
def string2index(stridx: str) -> Union[int, slice, str]:
    """Convert index string to either int or slice"""
    if ':' not in stridx:
        # may contain database accessor
        try:
            return int(stridx)
        except ValueError:
            return stridx
    i = [None if s == '' else int(s) for s in stridx.split(':')]
    return slice(*i)


def integer_to_slice_string(index):
    """
    Convert integer index to slice string.

    Args:
        index: Index to convert.

    Returns:
        slice: Converted slice.
    """
    if isinstance(index, slice):
        return index
    
    if is_convertible_to_integer(index):
        index=int(index)

    if isinstance(index, int):
        return string2index(f"{index}:{index+1}")
    elif index is None:
        return slice(None,None,None)
    elif isinstance(index, str):
        try:
            return string2index(index)
        except:
            raise ValueError("error creating slice from string {:s}".format(index))
    else:
        raise ValueError("`index` can be int, str, or slice, not {}".format(index))
    
#---------------------------------------#
def convert(what:Union[np.ndarray,float], family:str=None, _from:str="atomic_unit", _to:str="atomic_unit")->Union[np.ndarray,float]:
    """Convert a quantity from one unit of a specific family to another.
    Example: 
    arr = convert([1,3,4],'length','angstrom','atomic_unit')
    arr = convert([1,3,4],'energy','atomic_unit','millielectronvolt')"""
    # from ipi.utils.units import unit_to_internal, unit_to_user
    if family is not None:
        factor = unit_to_internal(family, _from, 1)
        factor *= unit_to_user(family, _to, 1)
        return what * factor
    else :
        return what
    
#---------------------------------------#
def get_offset(file:TextIOWrapper,
               Nmax:int=1000000,
               line_offset_old:list=None,
               index:slice=None):
    """
    Get line offsets in a file.

    Args:
        file (TextIOWrapper): File object.
        Nmax (int, optional): Maximum number of offsets. Defaults to 1000000.
        line_offset_old (list, optional): Old line offsets. Defaults to None.
        index (slice, optional): Slice index. Defaults to None.

    Returns:
        list: List of line offsets.
    """
    # Read in the file once and build a list of line offsets
    n = 0 
    line_offset = [None]*Nmax
    if line_offset_old is not None:
        line_offset[:len(line_offset_old)] = line_offset_old
        n = len(line_offset_old)
        del line_offset_old
    offset = 0
    restart = False
    if n == 0 : file.seek(0) # start from the beginning of the file
    for line in file: # cycle over all lines
        if n >= Nmax: # create a bigger list
            restart = True
            break
        if line.replace(" ","").startswith("#"): # check if the line contains a comment
            line_offset[n] = offset 
            n += 1
        if index is not None and index.stop is not None and n >= index.stop: # stop
            break
        offset += len(line)        
    if restart: return get_offset(file, Nmax=Nmax * 2, line_offset_old=line_offset, index=index)
    file.seek(0)
    return line_offset[:n]

#---------------------------------------#
def read_comments_xyz(file:TextIOWrapper,
                      index:slice=None):
    """
    Read comments from an XYZ file.

    Args:
        file (TextIOWrapper): File object.
        index (slice, optional): Slice index. Defaults to None.

    Returns:
        list: List of comments.
    """
    offset = get_offset(file=file,index=index)
    if index is not None:
        offset = offset[index]
    comments = [None]*len(offset)
    for n,l in enumerate(offset):
        file.seek(l)
        comments[n] = file.readline()
    return comments

#---------------------------------------#
class FakeList:
    """
    A fake list implementation.
    """
    def __init__(self, value: Any, length: int):
        self.value = value
        self.length = length

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, index: Union[int, slice]) -> Union[Any, List[Any]]:
        if isinstance(index, slice):
            start, stop, step = index.indices(self.length)
            return [self.value] * ((stop - start + step - 1) // step)
        elif 0 <= index < self.length:
            return self.value
        else:
            raise IndexError("FakeList index out of range")
        
#---------------------------------------#
def string2cell(string:Match[str])->np.ndarray:
    a, b, c = [float(x) for x in string.group(1).split()[:3]]
    alpha, beta, gamma = [float(x) * deg2rad for x in string.group(1).split()[3:6]]
    return abc2h(a, b, c, alpha, beta, gamma)

#---------------------------------------#
@dataclass
class Instruction:
    name:str
    bead:int
    filename:str
    format:str
    index:str
    pbc:bool
    fixed_cell:bool
        
    def get(self)->List[Atoms]:
        
        f = "extxyz" if self.format in ["i-pi","ipi"] else self.format
        
        units = {
            "positions" : "atomic_unit",
            "cell" : "atomic_unit"
        }
        # factor = {
        #     "positions" : np.nan,
        #     "cell" : np.nan
        # }
    
        with open(self.filename,"r") as ffile:
            traj = read(ffile,index=self.index,format=f)
            
            ffile.seek(0)
            
            if self.format in ["i-pi","ipi"]:
                
                # reading the comments
                comment = read_comments_xyz(ffile,slice(0,1,None))[0]
                
                # this could be improved
                units["positions"] = str(comment).split("positions{")[1].split("}")[0]
                units["cell"] = str(comment).split("cell{")[1].split("}")[0]
                
                factor = {
                    "positions" : convert(1,"length",_from=units["positions"],_to="angstrom"),
                            "cell" : convert(1,"length",_from=units["cell"],     _to="angstrom"),
                }

                if self.fixed_cell:
                    # little optimization if the cell is fixed
                    string:Match[str] = abcABC.search(comment)
                    value = string2cell(string)
                    cells:List[np.ndarray] = FakeList(value,len(traj))
                        
                else:
                    comments = read_comments_xyz(ffile,self.index)
                    strings:List[Match[str]] = [ abcABC.search(comment) for comment in comments ]
                    cells:List[np.ndarray] = [np.zeros(3,3)]*len(strings)
                    if len(comments) != len(traj):
                        raise ValueError("coding error: found comments different from atomic structures: {:d} comments != {:d} atoms (using index {})."\
                                            .format(len(comments),len(traj),self.index))
                    
                    # reading the cell parameters from the comments
                    for n,string in enumerate(strings):
                        cells[n] = string2cell(string)

                for atom, cell in zip(traj,cells):
                    atom.set_cell(cell.T * factor["cell"] if self.pbc else None)
                    atom.set_pbc(self.pbc)
                    atom.positions *= factor["positions"]
        
        return traj

    def __repr__(self):
        return f"Instruction(name={self.name}, bead={self.bead}, filename={self.filename}, format={self.format})"
    
    def __str__(self):
        return f"Instruction(name={self.name}, bead={self.bead}, filename={self.filename}, format={self.format})"
    
# class Instructions(List[Instruction]):
#     pass
#---------------------------------------#
def prepare_args():
    parser = argparse.ArgumentParser(description=description)
    argv = {"metavar":"\b"}
    parser.add_argument("-p", "--prefix"      , type=str     , required=False, **argv, help="prefix of the i-PI files (i-pi format will be assumed) (default: %(default)s)", default=None)
    parser.add_argument("-f", "--folder"      , type=str     , required=False, **argv, help="folder containing the i-PI files (default: %(default)s)", default=".")
    parser.add_argument("-n", "--index"       , type=itype   , required=False, **argv, help="index to be read from input file (default: %(default)s)", default=':')
    parser.add_argument("-pbc", "--pbc"       , type=str2bool, required=False, **argv, help="whether the output should be periodic (default: %(default)s)", default=True)
    parser.add_argument("-fc", "--fixed_cell" , type=str2bool, required=False, **argv, help="whether the cell is fixed along the trajectory (default: %(default)s)", default=False)
    parser.add_argument("-i", "--instructions", type=intype  , required=False, **argv, help="instructions as a list of tuples (name,filepath,format), e.g. [('positions', posfile, 'i-pi') , ('velocities', velfile, 'xyz'), ... ] (default: %(default)s)", default="automatic")
    parser.add_argument("-o", "--output"      , type=str     , required=True , **argv, help="output extxyz file")
    return parser.parse_args()

#---------------------------------------#
def main():
    
    #------------------#    
    args = prepare_args()    
    print(f"\n\t {description}")    
    print("\n\t Input arguments:")
    for k in args.__dict__.keys():
        print("\t {:>20s}:".format(k), getattr(args, k))
    print()
    
    #------------------#
    if args.prefix is None and args.instructions == "automatic":
        raise ValueError("Either a prefix or instructions must be given.")

    #------------------#
    if args.instructions != "automatic":
        raise ValueError("not implemented yet")

    #------------------#
    if args.prefix is not None:
        args.prefix = str(args.prefix)
        
    #------------------#  
    pattern = f"{args.folder}/{args.prefix}.*_*.xyz"
    print(f"\n\t Looking for files mathing the following pattern: '{pattern}'")
    files = glob.glob(pattern)
    N = len(files)
    print(f"\t {N} files found.")
    if N == 0:
        raise ValueError("No files found.")
    
    #------------------#
    pattern_extract = r'([^/]+)\.(.*?)_(.*?)\.(\w+)$'
    beads = set()
    names = set()
    for n,file in enumerate(files):
        matched = re.search(pattern_extract, os.path.basename(file))
        name = matched.group(2)    # Captures the first `*`
        bead = matched.group(3)   # Captures the second `*`
        names.add(name)
        beads.add(bead)  
    Nn = len(names)        
    Nb = len(beads)  
    print(f"\t {Nb} beads found.")    
    print(f"\t {Nn} arrays found.")    
    assert N == Nb*Nn, f"Found {N} files but {Nb} beads and {Nn} arrays."
        
    #------------------#
    print("\n\t Found the following files:")
    index = integer_to_slice_string(args.index)
    instructions = {}
    for n,file in enumerate(files):
        matched = re.search(pattern_extract, os.path.basename(file))
        if matched:
            prefix_value = matched.group(1)   # Captures the prefix (variable part)
            name = matched.group(2)    # Captures the first `*`
            bead = matched.group(3)   # Captures the second `*`
            extension = matched.group(4)      # Captures the file extension
            assert prefix_value == args.prefix, f"File {file} does not match the expected pattern"
            if extension != "xyz":
                warn(f"File {file} was expected to have extension 'xyz' but has extension '{extension}'")
            print(f"\t - bead={bead}, name={name}: {file}")
            if name not in instructions:
                instructions[name] = {}
            instructions[name][bead] = Instruction(name=name, bead=bead, filename=file, format="i-pi", index=index, fixed_cell=args.fixed_cell,pbc=args.pbc)
        else:
            raise ValueError(f"File {file} does not match the expected pattern")
        
    #------------------#
    print("\n\t Creating ASE trajectories (one for each bead)")
    trajs = [None]*Nb
    for n,bead in enumerate(beads):
        print(f"\t - bead={bead}")
        instr:Instruction = instructions["positions"][bead]
        trajs[n] = instr.get()
        pass
        # for name in names:
        #     for inst in instructions[name].values():
        #         if inst.bead == bead:
        
    pass

#---------------------------------------#
if __name__ == "__main__":
    main()
