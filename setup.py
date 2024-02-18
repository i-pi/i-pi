from pathlib import Path

from setuptools import setup, find_packages


setup(
    packages=[
        *find_packages(exclude=["ipi_tests*"]),
        "ipi._driver",
        "ipi._driver.pes",
    ],
    package_dir={
        "ipi._driver": "drivers/py",
    },
    scripts=[str(p) for p in Path("bin").iterdir()],
)
