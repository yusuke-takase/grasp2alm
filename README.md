# grasp2alm

This package supports the conversion from beam data calculated using [GRASP](https://www.ticra.com/software/grasp/) for CMB experiments to spherical harmonic coefficients ($a_{lm}$) based on the [HEALPix](https://healpix.sourceforge.io/) framework.
The code is designed based on [Beam](https://github.com/zonca/planck-levelS/tree/master/Beam), which is part of [LevelS](https://github.com/zonca/planck-levelS), the pipleline of the Planck experiment.

## Instllation

```
git clone https://github.com/yusuke-takase/grasp2alm
cd grasp2alm
pip install -e .
```

## Dependency

- `numpy`
- `healpy`
- `matplotlib`

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
- `to_map(nside, nstokes=3, outOftheta_val=0.0)`: Convert the beam to a map.
- `plot(stokes="I", color_resol=20, figsize=6, cmap="jet", return_fields=False)`: Plots the beam.

---

### Class: BeamMap

Represents a beam map.

#### Methods

- `to_alm(lmax, mmax)`: Converts the beam map to spherical harmonic coefficients, which is often called $b_{lm}$.

## Function Descriptions

### Function: grasp2alm

Convert a GRASP file to a spherical harmonic coefficients of `BeamMap`.

#### Arguments

- `filepath (str)`: Path to the GRASP file.
- `nside (int)`: nside used when generating a `BeamMap` from a GRASP file.
- `lmax (int)`: Maximum $l$ value for the spherical harmonic expansion.
- `mmax (int)`: Maximum $m$ value for the spherical harmonic expansion.
- `outOftheta_val (float)`: Value to assign to pixels outside the valid $\theta$ range in `BeamPolar`.
- `copol_axis (str, optional)`: Axis of the co-polarization component. Defaults to 'x'.
