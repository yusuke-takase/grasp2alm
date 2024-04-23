# grasp2alm

This package supports the conversion from beam data calculated using [GRASP](https://www.ticra.com/software/grasp/) for CMB experiments to spherical harmonic coefficients ($a_{lm}$) based on the [HEALPix](https://healpix.sourceforge.io/) framework.
The code is designed based on [Beam](https://github.com/zonca/planck-levelS/tree/master/Beam), which is part of [LevelS](https://github.com/zonca/planck-levelS), the pipleline of the Planck experiment.

## Instllation

```
git clone https://github.com/yusuke-takase/grasp2alm
cd grasp2alm
pip install -e .
```

## Class Descriptions

### Class: BeamGrid

Class to hold the data from a beam grid file of GRASP.

#### Methods

- `to_polar(self, copol_axis="x")`: Converts the beam grid to polar coordinates.
- `plot(self, pol='co', color_resol=20, figsize=6, cmap="jet", return_fields=False)`: Plots the beam grid.

---

### Class: BeamCut

Class to hold the data from a beam cut file of GRASP.

#### Methods

- `to_polar(self, copol_axis="x")`: Converts the beam to polar coordinates.
- `plot(self, pol='co', color_resol=20, figsize=6, cmap="jet",return_fields=False)`: Plots the beam.

---

### Class: BeamPolar

Represents beams of Stokes parameters.
Type to store Stokes parameters of a beam on a spherical polar ($\theta\phi$) grid.
This type is an input to the spherical harmonic transform.
Internally, this package uses the 'Ludwig 3' definition for thepolarisation basis with the co-polar (positive Q) direction aligned with the y-axis.

#### Methods

- `stokes_rotate()`: Rotates $Q$ and $U$ Stokes parameters from the co-cross basis to the polar basis.
  The $Q$ and $U$ Stokes parameters are usually represented in the
  co-cross basis, where the co-polar direction is aligned with the
  y-axis (consistent with Ludwig 3 convention). For the purposes of
  extracting the spherical harmonic coefficients, it is more useful
  to represent them in the polar basis. Unlike the [original Fortran
  method](https://github.com/zonca/planck-levelS/blob/master/Beam/beam_polar.f90#L260), this function operates on a deepcopy of the input `BeamPolar`,
  so it does not perform destructive actions.
- `to_map(nside, nstokes=3, outOftheta_val=0.0, interp_method="linear")`: Convert the beam to a map. Supported interp_methods are 'linear', 'nearest', 'slinear', 'cubic', 'quintic' and 'pchip'.
- `plot(stokes="I", color_resol=20, figsize=6, cmap="jet", return_fields=False)`: Plots the beam.

---

### Class: BeamMap

Represents a beam map.

#### Methods

- `to_alm(lmax=None,mmax=None,iter=3,pol=True,use_weights=False,datapath=None,gal_cut=0,use_pixel_weights=False,)`: Converts the beam map to spherical harmonic coefficients, which is often called $b_{lm}$.
- `to_alm_lsq(lmax,mmax,pol=True,tol=1e-10,maxiter=20)`
  Converts the beam map to spherical harmonic coefficients by using `healpy.map2alm_lsq`.

## Function Descriptions

### Function: grasp2alm

Convert a GRASP file to a spherical harmonic coefficients of `BeamMap` by using `healpy.map2alm`.

#### Arguments

- `filepath` (str): Path to the GRASP file.
- `nside` (int): Resolution parameter for the output beam map.
- `lmax` (int): The desired `lmax` parameters for the analysis.
- `mmax` (int): The desired `mmax` parameters for the analysis.
- `outOftheta_val` (float): Value to assign to pixels outside the valid `theta` range.
- `interp_method` (str): Interpolation method to use. Default is 'linear'. Supported are 'linear', 'nearest', 'slinear', 'cubic', 'quintic' and 'pchip'.
- `copol_axis` (str, optional): Axis of the co-polarization component. Defaults to 'x'.
- `iter` (int, scalar, optional): Number of iterations (default: 3).
- `pol` (bool, optional): If True, assumes input maps are TQU. Output will be TEB alm's.
  (input must be 1 or 3 maps)
  If False, apply spin 0 harmonic transform to each map.
  (input can be any number of maps)
  If there is only one input map, it has no effect. Default: True.
- `use_weights` (bool, scalar, optional): If True, use the ring weighting. Default: False.
- `datapath` (None or str, optional): If given, the directory where to find the pixel weights.
  See in the docstring above details on how to set it up.
- `gal_cut` (float [degrees]): Pixels at latitude in [-gal_cut;+gal_cut] are not taken into account.
- `use_pixel_weights` (bool, optional): If True, use pixel by pixel weighting, healpy will automatically download the weights, if needed.

### Function: grasp2alm_lsq

Convert a GRASP file to a spherical harmonic coefficients of `BeamMap` by using `healpy.map2alm_lsq`.

#### Arguments

- `filepath` (str): Path to the GRASP file.
- `nside` (int): Resolution parameter for the output beam map.
- `lmax` (int): The desired `lmax` parameters for the analysis.
- `mmax` (int): The desired `mmax` parameters for the analysis.
- `outOftheta_val` (float): Value to assign to pixels outside the valid `theta` range.
- `interp_method` (str): Interpolation method to use. Default is 'linear'. Supported are 'linear', 'nearest', 'slinear', 'cubic', 'quintic' and 'pchip'.
- `copol_axis` (str, optional): Axis of the co-polarization component. Defaults to 'x'.
- `pol` (bool, optional): If True, assumes input maps are TQU. Output will be TEB alm's.
  (input must be 1 or 3 maps)
  If False, apply spin 0 harmonic transform to each map.
  (input can be any number of maps)
  If there is only one input map, it has no effect. Default: True.
- `tol` (float): The desired accuracy for the result. Once this is reached, the iteration stops.
- `maxiter` (int): Maximum iteration count after which the minimization is stopped.
