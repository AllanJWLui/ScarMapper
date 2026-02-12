"""
Minimal setup.py — exists only to declare Cython extension modules.
All package metadata and dependencies live in pyproject.toml.

Delete scarmapper/setup.py; this is the only setup.py needed.
"""

from setuptools import setup, Extension
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(
        [
            Extension(
                "scarmapper.SlidingWindow",
                sources=["scarmapper/SlidingWindow.pyx"],
            ),
        ],
        language_level=3,
        annotate=False,
    ),
)
