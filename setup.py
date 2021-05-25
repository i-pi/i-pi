from pathlib import Path

from setuptools import setup, find_packages


setup(
    packages=[
        *find_packages(exclude=["ipi_tests*"]),
        "ipi._driver",
        "ipi._driver.pes",
        "ipi.tests",
        "ipi.tests.regression_tests",
    ],
    package_dir={
        "ipi._driver": "drivers/py",
        "ipi.tests": "ipi_tests",
    },
    package_data={"ipi.tests.regression_tests": ["tests/NVE/NVE_1/harmonic_python/*"]},
    scripts=[str(p) for p in Path("bin").iterdir()],
)
