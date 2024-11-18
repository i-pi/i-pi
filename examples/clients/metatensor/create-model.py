# Use https://github.com/metatensor/lj-test/ to define a basic LJ model,
# with and without a custom TorchScript extension
import metatensor_lj_test

# LJ parameters for Ni are taken from https://doi.org/10.1021/jp801931d
model = metatensor_lj_test.lennard_jones_model(
    atomic_type=28,
    cutoff=6.5,
    sigma=2.552,
    epsilon=5.65,
    length_unit="Angstrom",
    energy_unit="kcal/mol",
    with_extension=False,
)

model.save("nickel-lj.pt")


model = metatensor_lj_test.lennard_jones_model(
    atomic_type=28,
    cutoff=6.5,
    sigma=2.552,
    epsilon=5.65,
    length_unit="Angstrom",
    energy_unit="kcal/mol",
    with_extension=True,
)
model.save("nickel-lj-extensions.pt", collect_extensions="collected-extensions/")
