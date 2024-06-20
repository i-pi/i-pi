from ase.io import iread, write

position_filename = "../0_reference_pimd/simulation.xc.extxyz"
centroid_force_filename = "../0_reference_pimd/simulation.centroid_force.extxyz"
physical_force_filename = (
    "../0_reference_pimd/simulation.physical_force_on_centroid.extxyz"
)

atoms_list = []
for a_pos, a_cforce, a_pforce in zip(
    iread(position_filename),
    iread(centroid_force_filename),
    iread(physical_force_filename),
):

    atoms = a_pos.copy()
    atoms.set_positions(a_pos.arrays["x_centroid"])
    del a_pos.arrays["x_centroid"]

    atoms.arrays["centroid_force"] = a_cforce.arrays["forces_component_raw"]
    atoms.arrays["physical_force"] = a_pforce.arrays["forces_component_raw"]
    atoms.arrays["delta_force"] = (
        a_cforce.arrays["forces_component_raw"]
        - a_pforce.arrays["forces_component_raw"]
    )

    atoms_list.append(atoms.copy())

write("dataset.xyz", atoms_list)
