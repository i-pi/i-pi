from metatomic.testing import LennardJonesModel

# same parameters as the native plugin options in README.md
model = LennardJonesModel(
    atomic_type=18,
    cutoff=7.0,
    epsilon=0.0104,
    sigma=3.4,
    length_unit="angstrom",
    energy_unit="eV",
)
