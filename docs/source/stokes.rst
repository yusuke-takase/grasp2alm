Stokes conversion
=================

Conversion to Stokes parameter is performed following `Challinor et al. <https://arxiv.org/abs/astro-ph/0008228>`_ formalism.

Assuming a monochromatic system, the electric far field could be written as:

.. math::
    \vec{E} \propto \frac{1}{r}\Re{\{\vec{\mathcal{E}}\exp[i(kr-\omega t)]}\}

where :math:`\vec{\mathcal{E}}` is a complex, transverse vector function on the sphere, :math:`r` is the radial distance and :math:`\omega=ck` is the mean angular frequency of the radiation.

Writing the components of the field on the orthonormal basis vectors :math:`\{\hat\sigma_{\theta},\hat\sigma_{\phi}\}` of a spherical polar coordinate system as :math:`\mathcal{E}_{\theta}, \mathcal{E}_{\phi}` we can estimate the stokes parameter for the beam as:

.. math::
    I = \langle |\mathcal{E}_{\theta}|^2 + |\mathcal{E}_{\phi}|^2 \rangle \\
    Q = \langle |\mathcal{E}_{\theta}|^2 - |\mathcal{E}_{\phi}|^2 \rangle \\
    U = -2\Re\langle\mathcal{E}_{\theta}\mathcal{E}_{\phi}^*\rangle \\
    V = 2\Im\langle\mathcal{E}_{\theta}\mathcal{E}_{\phi}^*\rangle \\

The intensity I and circular polarization V are invariant under rotations of the polarization basis vectors, while Q and U transform like the components of a second-rank tensor.

Co- and Cross- polarized basis
------------------------------

The polar basis :math:`\{\hat\sigma_{\theta},\hat\sigma_{\phi}\}` is fixed relative to the sky and is singular at the north and south poles. For describing the beam, it is standard practice to use an alternative basis that is fixed relative to the detector and has its only singularity in the opposite direction to the main beam. 

We deﬁne a set of Cartesian basis vectors :math:`\{\hat\sigma_{x}',\hat\sigma_{y}', \hat\sigma_{z}'\}` which are ﬁxed relative to the detector. It is convenient to take :math:`\hat\sigma_{z}'` to be along the (nominal) main beam, while the polarization direction could be oriented along one of the other two axes.

Using this Cartesian frame we derive a set of polar basis vectors :math:`\{\hat\sigma_{\theta}',\hat\sigma_{\phi}'\}` on the sphere in the standard manner. The co- and cross-polar basis vectors are then derived as follows according to the `Ludwig-3 definition <https://ieeexplore.ieee.org/document/1140406>`_:

.. math::
    \hat\sigma_{co} = \sin\phi'\hat\sigma_{\theta}' + \cos\phi'\hat\sigma_{\phi}' \\
    \hat\sigma_{cross} = \cos\phi'\hat\sigma_{\theta}' - \sin\phi'\hat\sigma_{\phi}'

where :math:`\theta', \phi'` are spherical polar coordinates. 

To rotate from the co- and cross-polarized basis to the spherical polar basis we have to rotate through π/2 − φ′ in a right-handed sense about the inward normal to the sphere. Transforming the Stokes parameters for the beam from the co- and cross-polar basis to the spherical polar basis, we have:

.. math::
    I = \langle |\mathcal{E}_{co}|^2 + |\mathcal{E}_{cross}|^2 \rangle \\
    Q = -\langle |\mathcal{E}_{co}|^2 - |\mathcal{E}_{cross}|^2 \rangle \cos{2\phi'} + 2\Re\langle\mathcal{E}_{co}\mathcal{E}_{cross}^*\rangle \sin{2\phi'}\\
    U = -\langle |\mathcal{E}_{co}|^2 - |\mathcal{E}_{cross}|^2 \rangle \sin{2\phi'} - 2\Re\langle\mathcal{E}_{co}\mathcal{E}_{cross}^*\rangle \cos{2\phi'}\\
    V = -2\Im\langle\mathcal{E}_{\theta}\mathcal{E}_{\phi}^*\rangle \\

which reﬂect the spin-2 nature of the linear polarization.

After obtaining the Stokes parameters for the beam, the next step is to project them into a healpix map by interpolating the missing values. See :ref:`mapping:Conversion to healpix map`.