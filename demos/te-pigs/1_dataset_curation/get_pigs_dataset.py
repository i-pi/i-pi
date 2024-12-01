from ase.io import iread, write
from ipi import read_trajectory

position_filename = "../0_reference_pimd/simulation.xc.xyz"
centroid_force_filename = "../0_reference_pimd/simulation.centroid_force.extxyz"
physical_force_filename = "../0_reference_pimd/simulation.physical_force.extxyz"

atoms_list = []
for a_pos, a_cforce, a_pforce in zip(
    read_trajectory(position_filename),
    iread(centroid_force_filename),
    iread(physical_force_filename),
):
    atoms = a_pos.copy()

    atoms.arrays["centroid_force"] = a_cforce.arrays["f_centroid"]
    atoms.arrays["physical_force"] = a_pforce.arrays["forces_component_raw"]
    atoms.arrays["delta_force"] = (
        atoms.arrays["centroid_force"] - atoms.arrays["physical_force"]
    )

    atoms_list.append(atoms.copy())

write("dataset.xyz", atoms_list)
