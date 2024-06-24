# -*- encoding: utf-8 -*-

from dataclasses import dataclass, field
import numpy as np
import matplotlib.pyplot as plt
import warnings
from .beam_polar import BeamPolar

@dataclass
class BeamGrid:
    """Class to hold the data from a beam grid file of GRASP.

    Attributes:
        header (str): Record with identification text.
        ktype (int): Specifies type of file format.
        nset (int): Number of field sets or beams.
        klimit (int): Specification of limits in a 2D grid.
        icomp (int): Control parameter of field components.
        ncomp (int): Number of field components.
        igrid (int): Control parameter of field grid type.
        ix (int): Centre of set or beam No. i.
        iy (int): Centre of set or beam No. i.
        xs (float): Start x-coordinate of the grid.
        ys (float): Start y-coordinate of the grid.
        xe (float): End x-coordinate of the grid.
        ye (float): End y-coordinate of the grid.
        nx (int): Number of columns.
        ny (int): Number of rows.
        freq (float): Frequency.
        frequnit (str): Frequency unit.
        amp (np.ndarray): Array of complex amplitudes [theta, phi].

    Methods:
        __init__(self, filepath): Initialize the BeamGrid object.
        __post_init__(self): Read and parse the beam grid file.
        to_polar(self, copol_axis="x"): Convert the beam grid
            to polar coordinates.
        plot(self, pol='co', color_resol=20, figsize=6, cmap="inferno",
            return_fields=False): Plot the beam grid.

    """
    header: str = ""
    ktype: int = 0
    nset: int = 0
    klimit: int = 0
    icomp: int = 0
    ncomp: int = 0
    igrid: int = 0
    ix: int = 0
    iy: int = 0
    xs: float = 0.0
    ys: float = 0.0
    xe: float = 0.0
    ye: float = 0.0
    nx: int = 0
    ny: int = 0
    freq: float = 0.0
    frequnit: str = ""
    amp: np.ndarray = None


    def __init__(self, filepath):
        super().__init__()
        self.filepath = filepath
        self.filename = filepath.split("/")[-1]
        self.__post_init__()

    def __post_init__(self):
        if not self.filepath.endswith(".grd"):
            raise ValueError("Error in BeamGrid.__post_init__: The file is not a GRASP grid file.")
        with open(self.filepath, "r") as fi:
            while True:
                line = fi.readline()
                if line[0:4] == "++++":
                    break
                elif "FREQUENCIES [" in line:
                    self.frequnit = line.split("[")[1].split("]")[0]
                    self.freq = float(fi.readline().strip())
                else:
                    self.header += "\n" + line

            self.ktype = int(fi.readline())
            assert self.ktype == 1, "Unknown Grasp grid format, ktype != 1"

            line = fi.readline().split()
            self.nset  = int(line[0])
            self.icomp = int(line[1])
            self.ncomp = int(line[2])
            self.igrid = int(line[3])

            if self.nset > 1:
                warnings.warn("Warning: nset > 1, only reading first beam in file")

            line = fi.readline().split()
            self.ix = int(line[0])
            self.iy = int(line[1])
            i = 2
            while i <= self.nset:
                fi.readline()
                i += 1

            line = fi.readline().split()
            self.xs = float(line[0])
            self.ys = float(line[1])
            self.xe = float(line[2])
            self.ye = float(line[3])

            beam_solid_angle_rad = (np.cos(np.deg2rad(self.ys)) - np.cos(np.deg2rad(self.ye))) * (np.deg2rad(self.xe) - np.deg2rad(self.xs))
            if not np.isclose(beam_solid_angle_rad, 2.0*np.pi) and not np.isclose(beam_solid_angle_rad, 4.0*np.pi):
                warnings.warn(f"Warning: beam solid angle is not 2pi or 4pi because BeamGrid has xs={self.xs}, xe={self.xe}, ys={self.ys} and ye={self.ye}. The header should be checked.")

            line = fi.readline().split()
            self.nx = int(line[0])
            self.ny = int(line[1])
            self.klimit = int(line[2])
            self.amp = np.zeros((self.ncomp, self.nx, self.ny), dtype=complex)

            is_val, in_val = 0, 0
            for j in range(self.ny):
                if self.klimit == 1:
                    line = fi.readline().split()
                    in_val = int(line[1])
                    is_val = int(line[0])
                else:
                    is_val = 1
                    in_val = self.nx

                for i in range(is_val, in_val+1):
                    line = fi.readline().split()
                    self.amp[0, i-1, j] = float(line[0]) + float(line[1]) * 1j
                    self.amp[1, i-1, j] = float(line[2]) + float(line[3]) * 1j


    def to_polar(self, copol_axis="x"):
        """Converts beam in polar grid format into Stokes
        parameters on a polar grid. The value of copol
        specifies the alignment of the co-polar basis
        ('x' or 'y') of the input GRASP file.

        Args:
            copol_axis (str): The copolarization axis. Must be 'x' or 'y'.

        Returns:
            BeamPolar: The beam grid in polar coordinates.


        Raises:
            ValueError: If the beam is not in the supported GRASP grid format.

        """
        copol_axis = copol_axis.lower()
        assert self.ncomp == 2, "Error in BeamGrid.to_polar: beam is not in linear 'co' and 'cx' components"

        assert self.igrid == 7, "Error in BeamGrid.to_polar: beam is not on theta-phi grid"
        assert abs(self.xs) <= 1e-5, "Error in BeamGrid.to_polar: phi coordinates does not start at zero"
        assert abs(self.xe - self.xs - 360.0) <= 1e-5, "Error in BeamGrid.to_polar: phi range is not 360 degrees"
        assert copol_axis in ["x", "y"], "Error in BeamGrid.to_polar: copol_axis must be 'x' or 'y'"

        nphi = self.nx - 1
        ntheta = self.ny
        theta_rad_min = np.deg2rad(self.ys)
        theta_rad_max = np.deg2rad(self.ye)
        swaptheta = theta_rad_min > theta_rad_max
        if swaptheta:
            print("Warning: swapping theta direction")
            theta_rad_min = np.deg2rad(self.ye)
            theta_rad_max = np.deg2rad(self.ys)
        beam_polar = BeamPolar(nphi, ntheta, theta_rad_min, theta_rad_max, self.filename)

        if self.icomp == 3:
            if copol_axis == "x":
                sign = -1
            elif copol_axis == "y":

                sign = 1
            else:
                raise ValueError("Error in bm_grid2polar: unknown value for copol")
            c = self.amp[0, :-1, :]
            x = self.amp[1, :-1, :]
            modc2 = np.abs(c)**2
            modx2 = np.abs(x)**2
            acaxs = c * np.conj(x)
            beam_polar.stokes[0, :, :] = modc2 + modx2
            beam_polar.stokes[1, :, :] = sign * (modc2 - modx2)
            beam_polar.stokes[2, :, :] = sign * 2.0 * np.real(acaxs)
            beam_polar.stokes[3, :, :] = 2.0 * np.imag(acaxs)
            if swaptheta:
                beam_polar.stokes = beam_polar.stokes[:, :, ::-1]
        elif self.icomp == 9:
            c = self.amp[1, :-1, :]
            modc2 = np.abs(c)**2
            beam_polar.stokes[0, :, :] = modc2
            beam_polar.stokes[1, :, :] = 0.0
            beam_polar.stokes[2, :, :] = 0.0
            beam_polar.stokes[3, :, :] = 0.0
            if swaptheta:
                beam_polar.stokes = beam_polar.stokes[:, :, ::-1]
        else:
            raise ValueError("Error in grid2square: beam is not in supported grid sub-format")
        return beam_polar


    def plot(self, pol='co', color_resol=20, figsize=6, cmap="inferno", return_fields=False):
        """Plot the beam pattern.

        Args:
            pol (str): The polarization to plot. Must be either 'co' or 'cx'.
            color_resol (int): The number of color levels in the plot.
            figsize (int): The size of the figure.
            cmap (str): The colormap to use for the plot.
            return_fields (bool): Whether to return the x, y, and z values.

        Returns:
        If return_fields is False:
            None

        If return_fields is True:
            x (ndarray): The x values of the plot.
            y (ndarray): The y values of the plot.
            z (ndarray): The z values of the plot.

        """
        assert pol == 'co' or pol == 'cx', "Error in BeamGrid.plot: pol must be 'co' or 'cx'"
        dx = (self.xe - self.xs) / (self.nx - 1)
        dy = (self.ye - self.ys) / (self.ny - 1)

        xcen = dx * self.ix
        ycen = dy * self.iy
        xgrid = np.arange(xcen, self.nx) * dx
        ygrid = np.arange(ycen, self.ny) * dy
        phi_grid, theta_grid = np.deg2rad(np.meshgrid(xgrid, ygrid))
        x = np.rad2deg(theta_grid * np.cos(phi_grid))
        y = np.rad2deg(theta_grid * np.sin(phi_grid))
        if pol == 'co':
            z = np.real(self.amp[0])**2 + np.imag(self.amp[0])**2
            dBz = 10*np.log10(z)
        elif pol == 'cx':
            z = np.real(self.amp[1])**2 + np.imag(self.amp[1])**2
            dBz = 10*np.log10(z)

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
