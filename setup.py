"""
Setup script for SPUDSIM
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text() if readme_file.exists() else ""

setup(
    name="spudsim",
    version="2.2.0",
    description="USDA-ARS ACSL Potato Simulation Model",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="USDA-ARS Adaptive Cropping Systems Laboratory",
    author_email="david.fleisher@usda.gov",
    url="https://github.com/USDA-ARS-ACSL/SPUDSIM",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        # Core dependencies (minimal)
    ],
    extras_require={
        "hdf5": ["h5py>=3.0"],
        "netcdf": ["netCDF4>=1.5"],
        "all": ["h5py>=3.0", "netCDF4>=1.5"],
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    keywords="agriculture potato simulation crop-model",
)
