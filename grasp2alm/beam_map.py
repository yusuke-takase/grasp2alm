import numpy as np
import healpy as hp

class BeamMap:
    """Represents a beam map.

    Attributes:
        nside (int): The resolution parameter of the map.
        map (numpy.ndarray): The beam map data.

    Methods:
        to_alm(lmax, mmax): Converts the beam map to spherical harmonic coefficients.

    """

    def __init__(self, map):
        self.nside = hp.npix2nside(map.shape[1])
        self.map = map

    def to_alm(
        self,
        lmax=None,
        mmax=None,
        iter=3,
        pol=True,
        use_weights=False,
        datapath=None,
        gal_cut=0,
        use_pixel_weights=False,
        ):
        """Converts the beam map to spherical harmonic coefficients.

        Args:
            lmax (int): Maximum l value for the spherical harmonic expansion.
            mmax (int): Maximum m value for the spherical harmonic expansion.
            iter : int, scalar, optional
                Number of iteration (default: 3)
            pol : bool, optional
                If True, assumes input maps are TQU. Output will be TEB alm's.
                (input must be 1 or 3 maps)
                If False, apply spin 0 harmonic transform to each map.
                (input can be any number of maps)
                If there is only one input map, it has no effect. Default: True.
            use_weights : bool, scalar, optional
                If True, use the ring weighting. Default: False.
            datapath : None or str, optional
                If given, the directory where to find the pixel weights.
                See in the docstring above details on how to set it up.
            gal_cut : float [degrees]
                pixels at latitude in [-gal_cut;+gal_cut] are not taken into account
            use_pixel_weights: bool, optional
                If True, use pixel by pixel weighting, healpy will automatically download the weights, if needed

        Returns:
            numpy.ndarray: The spherical harmonic coefficients.

        Raises:
            AssertionError: If lmax is greater than 3*nside-1 or if mmax is greater than lmax.

        """

        assert self.map.shape[0] <=3 , "Error in BeamMap.to_alm: map has more than 3 Stokes parameters"

        alm = hp.map2alm(
            self.map,
            lmax=lmax,
            mmax=mmax,
            iter=iter,
            pol=pol,
            use_weights=use_weights,
            datapath=datapath,
            gal_cut=gal_cut,
            use_pixel_weights=use_pixel_weights
            )
        return alm

    def to_alm_lsq(
        self,
        lmax,
        mmax,
        pol=True,
        tol=1e-10,
        maxiter=20
        ):
        """Converts the beam map to spherical harmonic coefficients by using healpy.map2alm_lsq.

        Args:
            lmax, mmax : int
                The desired lmax and mmax parameters for the analysis
            pol : bool, optional
                If True, assumes input maps are TQU. Output will be TEB alm's.
            (input must be 1 or 3 maps)
                If False, apply spin 0 harmonic transform to each map.
            (input can be any number of maps)
                If there is only one input map, it has no effect. Default: True.
            tol : float
                The desired accuracy for the result. Once this is reached,
                the iteration stops.
            maxiter : int
                maximum iteration count after which the minimization is stopped

        Returns:
            numpy.ndarray: The spherical harmonic coefficients.

        Raises:
            AssertionError: If lmax is greater than 3*nside-1 or if mmax is greater than lmax.

        """

        assert self.map.shape[0] <=3 , "Error in BeamMap.to_alm: map has more than 3 Stokes parameters"

        alm = hp.map2alm_lsq(
            self.map,
            lmax,
            mmax,
            pol=pol,
            tol=tol,
            maxiter=maxiter
            )
        return alm
