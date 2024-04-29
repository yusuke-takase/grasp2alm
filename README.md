<p align="center">
  <h1>
  <img src="https://private-user-images.githubusercontent.com/83496454/326472496-956d5c54-1d93-4f3f-9d93-d69b3ae9be0b.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MTQ0MDA4NjAsIm5iZiI6MTcxNDQwMDU2MCwicGF0aCI6Ii84MzQ5NjQ1NC8zMjY0NzI0OTYtOTU2ZDVjNTQtMWQ5My00ZjNmLTlkOTMtZDY5YjNhZTliZTBiLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNDA0MjklMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjQwNDI5VDE0MjI0MFomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTAyMjUxMGYxYTY4NmNiOWI4NjhmY2VmZDJmODhlMjBmODE0Y2ZjYzlmY2MzN2E5MGE2ZTdmNjI0YWI3ODQ2ODImWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0JmFjdG9yX2lkPTAma2V5X2lkPTAmcmVwb19pZD0wIn0.yt8_PH_B17hRjsHbDkYsEQpguTQ3WsMcdP9AjV-IYoU" alt="Logo">
  </h1>
</p>

This package supports the conversion from beam data calculated using [GRASP](https://www.ticra.com/software/grasp/) for CMB experiments to spherical harmonic coefficients ($a_{lm}$) based on the [HEALPix](https://healpix.sourceforge.io/) framework.
The code is designed based on [Beam](https://github.com/zonca/planck-levelS/tree/master/Beam), which is part of [LevelS](https://github.com/zonca/planck-levelS), the pipleline of the Planck experiment.

## Instllation

```
pip install grasp2alm
```

Or you can install from source by:

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

- `to_alm(lmax=None,mmax=None,iter=3,pol=True,use_weights=False,datapath=None,gal_cut=0,use_pixel_weights=False)`: Converts the beam map to spherical harmonic coefficients, which is often called $b_{lm}$.
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
