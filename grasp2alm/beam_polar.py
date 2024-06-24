# -*- encoding: utf-8 -*-

from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
import copy
import healpy as hp
from .beam_map import BeamMap
from scipy.interpolate import RegularGridInterpolator

@dataclass
class BeamPolar:
    """Represents beams of Stokes parameters.

    Type to store Stokes parameters of a beam on a spherical
    polar(theta-phi) grid. This type is an input to the spherical
    harmonic transform. Internally, this package uses the
    'Ludwig 3' definition for thepolarisation basis with
    the co-polar (positive Q) direction aligned with the y-axis.

    Attributes:
        nphi (int): Number of phi values.
        ntheta (int): Number of theta values.
        theta_rad_min (float): Minimum theta value in radians.
        theta_rad_max (float): Maximum theta value in radians.
        stokes (ndarray): Array of shape (4, nphi, ntheta) representing the Stokes parameters (I,Q,U,V).
        filename (str): Name of the file.

    Methods:
        stokes_rotate(): Rotate the stokes beam.
        to_map(nside, nstokes=3, outOftheta_val=0.0): Convert the beam to a map.
        plot(stokes="I", color_resol=20, figsize=6, cmap="inferno", return_fields=False): Plot the beam.

    """

    def __init__(
        self,
        nphi: int,
        ntheta: int,
        theta_rad_min: float,
        theta_rad_max: float,
        filename: str
        ):
        self.nphi = nphi
        self.ntheta = ntheta
        self.theta_rad_min = theta_rad_min
        self.theta_rad_max = theta_rad_max
        self.stokes = np.zeros((4, nphi, ntheta), dtype=float)
        self.filename = filename

    def stokes_rotate(self):
        """Rotates Q and U Stokes parameters from the co-cross basis to the polar basis.
        The Q and U Stokes parameters are usually represented in the
        co-cross basis, where the co-polar direction is aligned with the
        y-axis (consistent with Ludwig 3 convention). For the purposes of
        extracting the spherical harmonic coefficients, it is more useful
        to represent them in the polar basis. Unlike the original LevelS's method, this function operates on a deepcopy of the input `BeamPolar`,
        so it does not perform destructive actions.

        Returns:
            BeamPolar: A new instance of BeamPolar with the rotated stokes beam.

        """
        beam_copy = copy.deepcopy(self)
        phi_step = 2 * np.pi / self.nphi
        theta_step = (self.theta_rad_max - self.theta_rad_min) / (self.ntheta - 1)

        theta_indices = np.arange(self.ntheta)
        theta_values = self.theta_rad_min + theta_indices * theta_step
        valid_theta_indices = theta_values != 0.0

        phi_indices = np.arange(self.nphi)
        phi_values = phi_indices * phi_step

        cos2phi = np.cos(2.0 * phi_values)
        sin2phi = np.sin(2.0 * phi_values)

        q = self.stokes[1, :, valid_theta_indices]
        u = self.stokes[2, :, valid_theta_indices]

        beam_copy.stokes[1, :, valid_theta_indices] =  q * cos2phi[None, :] + u * sin2phi[None, :]
        beam_copy.stokes[2, :, valid_theta_indices] = -q * sin2phi[None, :] + u * cos2phi[None, :]
        return beam_copy

    def to_map(self, nside, nstokes=3, outOftheta_val=hp.UNSEEN, interp_method="linear"):
        """Convert the BeamPolar to a BeamMap.

        Args:
            nside (int): The nside parameter for the HEALPix map.
            nstokes (int): Number of Stokes parameters.
            outOftheta_val (float): Value to fill outside the valid theta range.
            interp_method (str): Interpolation method to use. Default is 'linear'.
                Supported are 'linear', 'nearest', 'slinear', 'cubic', 'quintic' and 'pchip'.

        Returns:
            BeamMap: A new instance of BeamMap representing the beam map.

        """
        npix = hp.nside2npix(nside)
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        beam_polar = self.stokes_rotate()


        theta = theta[theta <= self.theta_rad_max]
        phi = phi[:len(theta)]

        beam_map = np.full((nstokes, npix), outOftheta_val, dtype=float)
        for s in range(nstokes):
            beam_map[s, :len(theta)] = beam_polar._get_interp_val(theta, phi, s, interp_method)
        return BeamMap(beam_map)

    def _get_interp_val(self, theta:np.ndarray, phi:np.ndarray, s:int, interp_method="linear"):
        """Calculate the value of the beam at a given theta, phi, and Stokes parameter.
        The value is interpolated from `BeamPolar` by a given theta and phi.

        Args:
            theta (float or array-like): The theta value(s) at which to evaluate the beam.
            phi (float or array-like): The phi value(s) at which to evaluate the beam.
            s (int): The Stokes parameter index.
            interp_method (str): Interpolation method to use. Default is 'linear'.
                Supported are 'linear', 'nearest', 'slinear', 'cubic', 'quintic' and 'pchip'.
        Returns:
            value (float or array-like): The value(s) of the beam at the given theta, phi, and Stokes parameter.

        """
        # Create a grid of theta and phi values
        theta_grid = np.linspace(self.theta_rad_min, self.theta_rad_max, self.ntheta)
        # To make a periodicity in phi, we add the first phi value to the end
        phi_grid = np.linspace(0.0, 2.0 * np.pi, self.nphi + 1)
        # To make a periodicity in stokes, we add the first stokes value to the end
        stokes_extended = np.concatenate([self.stokes[s], self.stokes[s][:1,:]], axis=0)

        # Create a 2D interpolator for the beam stokes values
        interpolator = RegularGridInterpolator(
            (phi_grid, theta_grid),
            stokes_extended,
            method=interp_method
            )

        # Use the interpolator to get the beam values at the given theta and phi
        value = interpolator(np.array([phi, theta]).T)

        return value


    def plot(self, stokes="I", color_resol=20, figsize=6, cmap="inferno", return_fields=False):
        """Plot the beam.

        Args:
            stokes (str): The Stokes parameter to plot. Can be "I", "Q", "U", or "V".
            color_resol (int): Number of color levels in the plot.
            figsize (int): Size of the figure.
            cmap (str): Colormap to use.
            return_fields (bool): Flag indicating whether to return the plot fields.

        Returns:
        If return_fields is False:
            None

        If return_fields is True:
            x (ndarray): The x values of the plot.
            y (ndarray): The y values of the plot.
            self.stokes[s] (ndarray): The z values of the plot.

        """
        if stokes.upper() =="I":
            s = 0
            label=r"$10\log_{10}(I)$ [dBi]"
        elif stokes.upper() == "Q":
            s = 1
            label=r"$10\log_{10}(|Q|)$ [dBi]"
        elif stokes.upper() == "U":
            s = 2
            label=r"$10\log_{10}(|U|)$ [dBi]"
        elif stokes.upper() == "V":
            s = 3
            label=r"$10\log_{10}(|V|)$ [dBi]"
        theta = np.linspace(self.theta_rad_min, self.theta_rad_max, self.ntheta)
        phi = np.linspace(0.0, 2.0 * np.pi, self.nphi)
        theta_grid, phi_grid = np.meshgrid(theta, phi)

        x = np.rad2deg(theta_grid * np.cos(phi_grid))
        y = np.rad2deg(theta_grid * np.sin(phi_grid))
        z = 10*np.log10(np.abs(self.stokes[s]))
        levels = np.linspace(np.min(z), np.max(z), color_resol)

        if not return_fields:
            plt.figure(figsize=(figsize,figsize))
            plt.axes().set_aspect("equal")
            plt.title(f"{self.filename} : ${stokes.upper()}$")
            plt.xlabel("Degrees")
            plt.ylabel("Degrees")
            plt.contourf(x, y, z, levels=levels, cmap=cmap, extend='both')
            plt.colorbar(orientation="vertical", label=label)
        else:
            return (x, y, self.stokes[s])


def _get_interp_val_from_polar_original(beam:BeamPolar, theta:np.ndarray, phi:np.ndarray, s:int):
    """Calculate the value of the beam at a given theta, phi, and Stokes parameter.
    The value is interpolated from `BeamPolar` by a given theta and phi.

    Args:
        beam (BeamPolar): The polar beam object.
        theta (float or array-like): The theta value(s) at which to evaluate the beam.
        phi (float or array-like): The phi value(s) at which to evaluate the beam.
        s (int): The Stokes parameter index.

    Returns:
        value (float or array-like): The value(s) of the beam at the given theta, phi, and Stokes parameter.

    """
    # TODO replace this function by scipy interpolation
    theta_step = (beam.theta_rad_max - beam.theta_rad_min) / (beam.ntheta - 1)
    phi_step = 2.0 * np.pi / beam.nphi
    theta = np.asarray(theta)
    phi = np.asarray(phi)

    ith1 = (theta / theta_step).astype(int)
    ith1 = np.maximum(0, np.minimum(beam.ntheta - 2, ith1))
    ith2 = ith1 + 1

    iph1 = (phi / phi_step).astype(int)
    iph1 = np.maximum(0, np.minimum(beam.nphi - 1, iph1))
    iph2 = iph1 + 1
    iph2[iph2 >= beam.nphi] = 0

    th1 = beam.theta_rad_min + ith1 * theta_step
    wth = 1.0 - (theta - th1) / theta_step

    ph1 = iph1 * phi_step
    wph = 1.0 - (phi - ph1) / phi_step

    value = wth * (wph * beam.stokes[s, iph1, ith1] + (1.0 - wph) * beam.stokes[s, iph2, ith1]) + \
            (1.0 - wth) * (wph * beam.stokes[s, iph1, ith2] + (1.0 - wph) * beam.stokes[s, iph2, ith2])
    return value
