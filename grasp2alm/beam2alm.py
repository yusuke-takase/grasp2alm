import numpy as np
import healpy as hp
from .beam_grid import BeamGrid
from .beam_cut import BeamCut

def grasp2alm(
    filepath,
    nside,
    outOftheta_val=hp.UNSEEN,
    interp_method="linear",
    copol_axis='x',
    lmax=None,
    mmax=None,
    iter=3,
    pol=True,
    use_weights=False,
    datapath=None,
    gal_cut=0,
    use_pixel_weights=False,
    ):
    """Convert a GRASP file to a spherical harmonic coefficients of beam map.

    Args:
        filepath (str): Path to the GRASP file.
        nside (int): Resolution parameter for the output beam map.
        lmax (int): The desired lmax parameters for the analysis.
        mmax (int): The desired mmax parameters for the analysis.
        outOftheta_val (float): Value to assign to pixels outside
            the valid theta range.
        interp_method (str): Interpolation method to use. Default is 'linear'.
                Supported are 'linear', 'nearest', 'slinear', 'cubic', 'quintic' and 'pchip'.
        copol_axis (str, optional): Axis of the co-polarization
            component. Defaults to 'x'.
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
        alm (numpy.ndarray): Spherical harmonic coefficients of the beam map.

    Raises:
        ValueError: If the file format is unknown.

    """
    if filepath.endswith(".grd"):
        beam = BeamGrid(filepath)
    elif filepath.endswith(".cut"):
        beam = BeamCut(filepath)
    else:
        raise ValueError("Error in grasp2alm: unknown file format")
    beam_polar = beam.to_polar(copol_axis)
    beam_map = beam_polar.to_map(
        nside,
        outOftheta_val=outOftheta_val,
        interp_method=interp_method
        )
    alm = beam_map.to_alm(
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

def grasp2alm_lsq(
    filepath,
    nside,
    lmax,
    mmax,
    outOftheta_val=hp.UNSEEN,
    interp_method="linear",
    copol_axis='x',
    pol=True,
    tol=1e-10,
    maxiter=20
    ):
    """Convert a GRASP file to a spherical harmonic coefficients of beam map by using healpy.map2alm_lsq.

    Args:
        filepath (str): Path to the GRASP file.
        nside (int): Resolution parameter for the output beam map.
        lmax (int): The desired lmax parameters for the analysis.
        mmax (int): The desired mmax parameters for the analysis.
        outOftheta_val (float): Value to assign to pixels outside
            the valid theta range.
        interp_method (str): Interpolation method to use. Default is 'linear'.
                Supported are 'linear', 'nearest', 'slinear', 'cubic', 'quintic' and 'pchip'.
        copol_axis (str, optional): Axis of the co-polarization
            component. Defaults to 'x'.
        pol : bool, optional
            If True, assumes input maps are TQU. Output will be TEB alm's.
            (input must be 1 or 3 maps)
            If False, apply spin 0 harmonic transform to each map.
            (input can be any number of maps)
            If there is only one input map, it has no effect. Default: True.
        tol : float
            The desired accuracy for the result. Once this is reached, the iteration
            stops.
        maxiter : int
            maximum iteration count after which the minimization is stopped

    Returns:
        alm (numpy.ndarray): Spherical harmonic coefficients of the beam map.

    Raises:
        ValueError: If the file format is unknown.

    """
    if filepath.endswith(".grd"):
        beam = BeamGrid(filepath)
    elif filepath.endswith(".cut"):
        beam = BeamCut(filepath)
    else:
        raise ValueError("Error in grasp2alm: unknown file format")
    beam_polar = beam.to_polar(copol_axis)
    beam_map = beam_polar.to_map(
        nside,
        outOftheta_val=outOftheta_val,
        interp_method=interp_method
        )
    alm = beam_map.to_alm_lsq(
        lmax,
        mmax,
        pol=pol,
        tol=tol,
        maxiter=maxiter
        )
    return alm
