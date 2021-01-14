from ase.io import read, write
from ase.geometry.analysis import Analysis
import numpy as np
import os



def set_pbc(pos, cell):
    for single_pos in pos:
        single_pos.set_cell(cell)
        single_pos.set_pbc(True)

def take_frame_rdf(pos, r, nbin, frames, elements):
    tmp_info = Analysis(pos)
    tmp_rdf_list = tmp_info.get_rdf(r, nbin, imageIdx=slice(0, frames, 1), elements=elements)
    gr_list = []
    for s_gr in tmp_rdf_list:
        gr_list.append(s_gr)
    gr_list = np.array(gr_list)
    gr_list = gr_list.T
    return gr_list

def xyz2rdflist(in_name, out_name, cell, r, nbin, elements, frames = None):
    pos = read(in_name, index= ":")
    set_pbc(pos, cell)
    if frames == None:
        frames = len(pos)
    rdflist = take_frame_rdf(pos, r, nbin, frames, elements)
    r_list = np.arange(0, r, r/nbin)
    tmp = np.concatenate((np.array([r_list]).T, rdflist), axis=1)
    np.savetxt(out_name, tmp.T[1:])
    np.savetxt("r_value.dat", tmp[0])

## user input area
in_file_name = ["simulation.pos_0.xyz"]
out_file_name = ["rdf_list_O-O.dat"]

cell = [35.23300, 35.23300, 35.23300, 90, 90, 90]
r = 7
nbin = 200
elements = ["O", "O"]
## user input area

for i in range(1):
    xyz2rdflist(in_file_name[i], out_file_name[i], cell, r, nbin, elements)
    print("write", out_file_name[i], "finished")
