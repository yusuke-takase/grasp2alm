# -*- encoding: utf-8 -*-

from dataclasses import dataclass, field
import numpy as np
import matplotlib.pyplot as plt
from typing import Union
from .beam_polar import BeamPolar

@dataclass
class BeamCut:
    """
    Class to hold the data from a beam cut file of GRASP.

    Args:
        header (str): Record with identification text.
        vini (float): Initial value.
        vinc (float): Increment.
        vnum (int): Number of values in cut.
        c (np.ndarray): Constant.
        icomp (int): Polarization control parameter.
        icut (int): Control parameter of cut.
        ncomp (int): Number of field components.
        ncut (int): Number of cuts.
        amp (np.ndarray): Amplitude.
    """
    header: str = ""
    vini: float = 0.0
    vinc: float = 0.0
    vnum: int = 0
    c = np.array([])

    icomp: int = 0
    icut: int = 0
    ncomp: int = 0
    ncut: int = 0
    amp: np.ndarray = None

    def __init__(self, filepath):
        """
        Initializes a BeamCut object.
        """
        super().__init__()
        self.filepath = filepath
        self.filename = filepath.split("/")[-1]
        self.__post_init__()

    def __post_init__(self):
        """
        Performs post-initialization tasks.
        """
        if not self.filepath.endswith(".cut"):
            raise ValueError("Error in BeamCut.__post_init__: The file is not a GRASP cut file.")
        with open(self.filepath, "r") as fi:
            self.header = fi.readline().strip()
            while True:
                line = fi.readline()
                if not line:
                    break
                data = line.split()
                if len(data) == 7:
                    self.vini, self.vinc, self.vnum, c, self.icomp, self.icut, self.ncomp = map(float, data)
                    self.vnum, self.icomp, self.icut, self.ncomp = map(int, (self.vnum, self.icomp, self.icut, self.ncomp))
                    self.c = np.append(self.c, c)
                    self.ncut += 1
                if self.ncomp > 2:
                    raise ValueError("Three field components present. Beam package can only handle two field components.")
                if self.vnum % 2 == 0:
                    raise ValueError("The number of pixels in a cut must be odd.")

            self.amp = np.zeros((self.ncomp, self.vnum, self.ncut), dtype=complex)
            fi.seek(0)
            cnt = 0
            while True:
                line = fi.readline()
                if not line:
                    break
                data = line.split()
                if len(data) == 7:
                    self.vini, self.vinc, self.vnum, _, self.icomp, self.icut, self.ncomp = map(float, data)
                    self.vnum, self.icomp, self.icut, self.ncomp = map(int, (self.vnum, self.icomp, self.icut, self.ncomp))
                    for i in range(self.vnum):
                        line = fi.readline()
                        data = line.split()
                        tmp1, tmp2, tmp3, tmp4 = map(float, data)
                        if any(np.isnan([tmp1, tmp2, tmp3, tmp4])):
                            raise ValueError("Encountered a NaN value in Amplitude. Please check your input.")
                        self.amp[0, i, cnt] = complex(tmp1, tmp2)
                        self.amp[1, i, cnt] = complex(tmp3, tmp4)
                    cnt += 1

    def to_polar(self, copol_axis="x"):
        """Converts beam in "cut" format to Stokes parameters
        on a polar grid.  Assumes that cuts are evenly spaced
        in theta. The value of copol specifies the alignment
        of the co-polar basis ('x' or 'y') of the input GRASP file.

        Args:
            copol_axis (str): The axis of copolarization. Must be either 'x' or 'y'.

        Returns:
            BeamPolar: The beam in polar coordinates.

        Raises:
            ValueError: If the beam is not in the expected format.

        """
        copol_axis = copol_axis.lower()

        if self.icomp != 3:
            raise ValueError("Error in BeamCut.to_polar: beam is not in linear 'co' and 'cx' components")
        if self.icut != 1:
            raise ValueError("Error in BeamCut.to_polar: beam is not in phi cuts")
        if self.ncomp != 2:
            raise ValueError("Error in BeamCut.to_polar: beam has the wrong number of components")
        if copol_axis not in ["x", "y"]:
            raise ValueError("Error in BeamCut.to_polar: copol_axis must be 'x' or 'y'")

        nphi = int(2 * self.ncut)
        ntheta = int(self.vnum // 2)
        theta_rad_min = 0.0
        theta_rad_max = np.deg2rad(np.abs(self.vini))
        beam_polar = BeamPolar(nphi, ntheta, theta_rad_min, theta_rad_max, self.filename)
        amp_tmp = np.zeros((2, nphi, ntheta), dtype=complex)

        for i in range(self.ncut):
            amp_tmp[:, i, :] = self.amp[:, ntheta:self.vnum-1, i]
            amp_tmp[:, self.ncut + i, :] = self.amp[:, ntheta-1::-1, i]

        if copol_axis == "x":
            sign = -1
        elif copol_axis == "y":
            sign = 1

        c = amp_tmp[0, :, :]
        x = amp_tmp[1, :, :]

        modc2 = np.abs(c)**2
        modx2 = np.abs(x)**2
        acaxs = c * np.conj(x)

        beam_polar.stokes[0, :, :] = modc2 + modx2
        beam_polar.stokes[1, :, :] = sign * (modc2 - modx2)
        beam_polar.stokes[2, :, :] = sign * 2.0 * np.real(acaxs)
        beam_polar.stokes[3, :, :] = 2.0 * np.imag(acaxs)
        return beam_polar

    def plot(
        self, 
        pol='co', 
        color_resol=20, 
        figsize=6, 
        cmap="jet", 
        return_fields=False
        ):
        """Plot the beam pattern.

        Args:
            pol (str): The polarization to plot. Must be either 'co' or 'cx'.
            color_resol (int): The number of color levels in the plot.
            figsize (int): The size of the figure.
            cmap (str): The colormap to use for the plot.
            return_fields (bool): Whether to return the x, y, and z values.

        Returns
        -------
            None: if return_fields is False (default)
            (ndarray,ndarray,ndarray): if return_fields is True returns x,y,z values of the plot.
        """
        if not (pol == 'co' or pol == 'cx'):
            raise ValueError("Error in BeamCut.plot: pol must be 'co' or 'cx'")

        theta = np.deg2rad(np.linspace(self.vini, self.vini + self.vinc * (self.vnum - 1), self.vnum))
        phi = np.deg2rad(np.linspace(0.0, 180.0, self.ncut))

        theta_grid, phi_grid = np.meshgrid(theta, phi)
        if pol == 'co':
            z = np.real(self.amp[0])**2 + np.imag(self.amp[0])**2
            dBz = 10*np.log10(z)
        elif pol == 'cx':
            z = np.real(self.amp[1])**2 + np.imag(self.amp[1])**2
            dBz = 10*np.log10(z)
        x = np.rad2deg(theta_grid * np.cos(phi_grid))
        y = np.rad2deg(theta_grid * np.sin(phi_grid))
        levels = np.linspace(np.min(dBz), np.max(dBz), color_resol)

        if not return_fields:
            plt.figure(figsize=(figsize,figsize))
            plt.axes().set_aspect("equal")
            plt.title(f"{self.filename} : {pol}")
            plt.xlabel("Degrees")
            plt.ylabel("Degrees")
            plt.contourf(x, y, dBz.T, levels=levels, cmap=cmap, extend='both')
            plt.colorbar(orientation="vertical", label="dBi")
        else:
            return (x, y, z.T)
