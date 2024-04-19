# levels_beam

This repository contains Python implementation of [the beam function for Planck LevelS](https://github.com/zonca/planck-levelS/tree/master/Beam).

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
Type to store Stokes parameters of a beam on a spherical polar (theta-phi) grid.
This type is an input to the spherical harmonic transform.
Internally, this package uses the 'Ludwig 3' definition for thepolarisation basis with the co-polar (positive Q) direction aligned with the y-axis.

#### Methods

- `stokes_rotate(healpy_pol_convention=True)`: Rotates Q and U Stokes parameters from the co-cross basis to the polar basis.
  The Q and U Stokes parameters are usually represented in the
  co-cross basis, where the co-polar direction is aligned with the
  y-axis (consistent with Ludwig 3 convention). For the purposes of
  extracting the spherical harmonic coefficients, it is more useful
  to represent them in the polar basis. This routine should only be
  called just before the spherical transform routines.
- `to_map(nside, healpy_pol_convention=True, nstokes=3, outOftheta_val=0.0)`: Convert the beam to a map.
- `plot(stokes="I", color_resol=20, figsize=6, cmap="jet", return_fields=False)`: Plots the beam.

---

### Class: BeamMap

Represents a beam map.

#### Methods

- `to_blm(lmax, mmax)`: Converts the beam map to spherical harmonic coefficients.

## Function Descriptions

### Function: grasp2blm

Convert a GRASP file to a spherical harmonic coefficients of beam map.

#### Arguments

- `filepath (str)`: Path to the GRASP file.
- `nside (int)`: Resolution parameter for the output beam map.
- `lmax (int)`: Maximum $l$ value for the spherical harmonic expansion.
- `mmax (int)`: Maximum $m$ value for the spherical harmonic expansion.
- `outOftheta_val (float)`: Value to assign to pixels outside the valid theta range.
- `copol_axis (str, optional)`: Axis of the co-polarization component. Defaults to 'x'.
- `healpy_pol_convention (bool, optional)`: Whether to use the Healpix polarization convention. Defaults to `True`.

---

### Function: \_get_beam_polar_value

#### Arguments

- `beam (BeamPolar)`: The polar beam object.
- `theta (float or array-like)`: The theta value(s) at which to evaluate the beam.
- `phi (float or array-like)`: The phi value(s) at which to evaluate the beam.
- `s (int)`: The Stokes parameter index.

#### Returns

- `value (float or array-like)`: The value(s) of the beam at the given theta, phi, and Stokes parameter.
