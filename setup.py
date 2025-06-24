from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize
import numpy as np
import glob
import os

# This is a stupid way of doing this VVV
extensions = cythonize(
    [
        Extension(
            "*",
            sources=["src/UM_MRF/**/*.py"],
            include_dirs=["src"],
        )
    ],
    compiler_directives={"language_level": "3"}
)

setup(
    name="UM_MRF",
    ext_modules=extensions,
    packages=["UM_MRF", "UM_MRF.sim_blocks"],
    package_dir={"": "src"},
    version="1.0",
    author="Christopher Louly",
    author_email="clouly@umich.edu",
    cmdclass={"build_ext": build_ext},
    zip_safe=False
)
