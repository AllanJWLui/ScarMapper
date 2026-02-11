"""
Setup file to Cythonize Alignment Processing using "python3 setup.py build_ext --inplace"
"""
import os
import pathlib
from distutils.core import setup
from Cython.Build import cythonize

# Only run setup when this script is executed directly
if __name__ == "__main__":
    setup(
        name="scarmapper",
        author='Dennis Simpson',
        author_email='dennis@email.unc.edu',
        ext_modules=cythonize("scarmapper/SlidingWindow.pyx", annotate=False)  # Note: path includes scarmapper/
    )
