.. grasp2alm documentation master file, created by
   sphinx-quickstart on Mon May 27 15:48:41 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

User's Manual of grasp2alm
==========================

This package supports the conversion from beam data calculated using `GRASP <https://www.ticra.com/software/grasp/>`_ for CMB experiments to spherical harmonic coefficients :math:`a_{\ell m}` based on the `HEALPix <https://healpix.sourceforge.io/>`_ framework.
The code is designed based on `Beam <https://github.com/zonca/planck-levelS/tree/master/Beam/>`_, which is part of `LevelS <https://github.com/zonca/planck-levelS>`_, the pipeline of the Planck experiment.

GitHub Repository
-----------------

You can find the source code and contribute to the project on `GitHub <https://github.com/yusuke-takase/grasp2alm>`_.

Contents
--------

.. toctree::
   :maxdepth: 1

   installation
   tutorial
   cut
   grid
   stokes
   mapping
   gauss_beam
   reference

Change Logs
-----------

Review the changes in each release in the `CHANGELOG on Github <https://github.com/yusuke-takase/grasp2alm/blob/master/CHANGELOG.md>`_.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
