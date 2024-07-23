Conversion to HEALPix map
=========================

To achieve the conversion into spherical harmonics it is necessary to project the beam into a `HELAPix map <https://en.wikipedia.org/wiki/HEALPix>`_, this allows grasp2alm to use built-in functions such as the `spherical harmonics transform <https://healpy.readthedocs.io/en/latest/healpy_spht.html>`_.

HELAPix map produces a subdivision of a spherical surface in which each pixel covers the same surface area as every other pixel. The resolution is controlled by the parameter :math:`N_{\mathrm{side}}` related to the number of pixels by the simple equation :math:`N_{\mathrm{pixel}}=12N_{\mathrm{side}}^2`.

The discretization performed by Grasp is different from the one used in a HEALPix map, to avoid gaps in the map it is necessary to make an interpolation of the beam at a given direction.

The interpolation is performed using `scipy.interpolate.RegularGridInterpolator <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html>`_; starting from the regular grid :math:`\theta, \phi` (or :math:`u,v`) returned by the reading and rotating step, the map is created interpolating the values returned by using `healpy.pix2ang <https://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.pix2ang.html>`_.

The interpolation could be performed by choosing between the following methods: `linear`, `nearest`, `slinear`, `cubic`, `quintic` and `pchip`. Moreover, it is possible to choose the value of the pixels that are over the :math:`\theta, \phi` boundaries, default is `healpy.pixelfunc.UNSEEN <https://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.UNSEEN.html>`_.