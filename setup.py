from pathlib import Path
from setuptools import setup, find_packages

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="pynoa",
    version='0.1.0',
    author="Bagaskara Primastya Putra",
    author_email="primastyaputra@gmail.com",
    description="Python Nonlinear Observability Analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BagaskaraPutra/PyNOA",
    packages=find_packages(),
    install_requires=[
        "ipykernel",
        "datetime",
        "IPython",
        "sympy",
        "dill",
    ],
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Software Development",
        "Topic :: Software Development :: Libraries",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Programming Language :: Python :: 3",
    ],
    keywords="observability analysis, control system",
)