from setuptools import setup
from levels_beam.version import __author__, __original_url__, __version__

setup(
    name = 'levels_beam',
    version = __version__,
    description = 'A package that converts Planck-LevelS beam-related code to python, supporting the computation of Stokes parameter beams from the GRASP beam format, projection onto the HEALPix sphere, conversion to expansion coefficients for spherical harmonics, etc.',
    author = __author__,
    author_email = 'takase_y@s.okayama-u.ac.jp',
    packages = ['levels_beam'],
    install_requires=[
        'numpy',
        'healpy',
        'matplotlib',
    ],
    license='MIT',
)
