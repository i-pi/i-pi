from pathlib import Path

from setuptools import setup, find_packages

# exclude i-pi-driver that is a compiled fortran file when installing from a local dir
scripts = [str(p) for p in Path("bin").iterdir() if p.name != "i-pi-driver"]

setup(
    packages=[*find_packages(exclude=["ipi_tests*"])],
    package_data={
        "ipi.utils": ["lebedev_grids.npy"],  # Include the npy file here
    },
    scripts=scripts,
)
