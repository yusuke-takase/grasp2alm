from dataclasses import dataclass
import numpy as np
import warnings
import copy
import healpy as hp
import matplotlib.pyplot as plt

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
        plot(self, pol='co', color_resol=20, figsize=6, cmap="jet",
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
    amp: np.ndarray = np.array([])


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
                print("Warning: nset > 1, only reading first beam in file")

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


    def plot(self, pol='co', color_resol=20, figsize=6, cmap="jet", return_fields=False):
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
@dataclass
class BeamCut:
    """Class to hold the data from a beam cut file of GRASP.

    Attributes:
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

    Methods:
        __init__(self, filepath): Initializes a BeamCut object.
        __post_init__(self): Performs post-initialization tasks.
        to_polar(self, copol_axis="x"): Converts the beam
            to polar coordinates.
        plot(self, pol='co', color_resol=20, figsize=6, cmap="jet",
            return_fields=False): Plots the beam.

    """
    header: str = ""   # TEXT - Record with identification text
    vini: float = 0.0  # Initial value
    vinc: float = 0.0  # Increment
    vnum: int = 0      # Number of values in cut
    c: np.ndarray = np.array([]) # Constant
    icomp: int = 0     # Polarization control parameter.
    icut: int = 0      # Control parameter of cut.
    ncomp: int = 0     # Number of field components
    ncut: int = 0      # Number of cuts
    amp: np.ndarray = np.array([])

    def __init__(self, filepath):
        super().__init__()
        self.filepath = filepath
        self.filename = filepath.split("/")[-1]
        self.__post_init__()

    def __post_init__(self):
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
                    self.vini, self.vinc, self.vnum, c, self.icomp, self.icut, self.ncomp = map(float, data)
                    self.vnum, self.icomp, self.icut, self.ncomp = map(int, (self.vnum, self.icomp, self.icut, self.ncomp))
                    for i in range(self.vnum):
                        line = fi.readline()
                        data = line.split()
                        tmp1, tmp2, tmp3, tmp4 = map(float, data)
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
        # Check metadata of input beam.
        if self.icomp != 3:
            raise ValueError("Error in BeamCut.to_polar: beam is not in linear 'co' and 'cx' components")
        if self.icut != 1:
            raise ValueError("Error in BeamCut.to_polar: beam is not in phi cuts")
        if self.ncomp != 2:
            raise ValueError("Error in BeamCut.to_polar: beam has the wrong number of components")
        assert copol_axis in ["x", "y"], "Error in BeamCut.to_polar: copol_axis must be 'x' or 'y'"

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

    def plot(self, pol='co', color_resol=20, figsize=6, cmap="jet", return_fields=False):
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
        assert pol == 'co' or pol == 'cx', "Error in BeamCut.plot: pol must be 'co' or 'cx'"

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
        stokes (ndarray): Array of shape (4, nphi, ntheta) representing the Stokes parameters.
        filename (str): Name of the file.

    Methods:
        stokes_rotate(healpy_pol_convention=True): Rotate the polarization beam.
        to_map(nside, healpy_pol_convention=True, nstokes=3, outOftheta_val=0.0): Convert the beam to a map.
        plot(stokes="I", color_resol=20, figsize=6, cmap="jet", return_fields=False): Plot the beam.

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

    def stokes_rotate(self, healpy_pol_convention=True):
        """Rotates Q and U Stokes parameters from the co-cross basis to the polar basis.

        The Q and U Stokes parameters are usually represented in the
        co-cross basis, where the co-polar direction is aligned with the
        y-axis (consistent with Ludwig 3 convention).  For the purposes of
        extracting the spherical harmonic coefficients, it is more useful
        to represent them in the polar basis.  This routine should only be
        called just before the spherical transform routines.

        Args:
            healpy_pol_convention (bool): Flag indicating whether to use the healpy polarization convention.

        Returns:
            BeamPolar: A new instance of BeamPolar with the rotated polarization beam.

        """
        beam_copy = copy.deepcopy(self)
        phi_step = 2 * np.pi / self.nphi
        theta_step = (self.theta_rad_max - self.theta_rad_min) / (self.ntheta - 1)

        if healpy_pol_convention == False:
            warnings.warn("Warning: Not healpy convention. Flip positive and negative value of polarization beam. See https://healpy.readthedocs.io/en/latest/blm_gauss_plot.html")

        theta_indices = np.arange(self.ntheta)
        theta_values = self.theta_rad_min + theta_indices * theta_step
        valid_theta_indices = theta_values != 0.0

        phi_indices = np.arange(self.nphi)
        phi_values = phi_indices * phi_step

        cos2phi = np.cos(2.0 * phi_values)
        sin2phi = np.sin(2.0 * phi_values)

        q = self.stokes[1, :, valid_theta_indices]
        u = self.stokes[2, :, valid_theta_indices]
        if healpy_pol_convention == False:
            q = -q
            u = -u

        beam_copy.stokes[1, :, valid_theta_indices] =  q * cos2phi[None, :] + u * sin2phi[None, :]
        beam_copy.stokes[2, :, valid_theta_indices] = -q * sin2phi[None, :] + u * cos2phi[None, :]
        return beam_copy

    def to_map(self, nside, healpy_pol_convention=True, nstokes=3, outOftheta_val=0.0):
        """Convert the BeamPolar to a BeamMap.

        Args:
            nside (int): The nside parameter for the HEALPix map.
            healpy_pol_convention (bool): Flag indicating whether to use the healpy polarization convention.
            nstokes (int): Number of Stokes parameters.
            outOftheta_val (float): Value to fill outside the valid theta range.

        Returns:
            BeamMap: A new instance of BeamMap representing the beam map.

        """
        npix = hp.nside2npix(nside)
        beampolar = self.stokes_rotate(healpy_pol_convention=healpy_pol_convention)
        theta, phi = hp.pix2ang(nside, np.arange(npix))

        theta = theta[theta <= self.theta_rad_max]
        phi = phi[:len(theta)]

        beam_map = np.full((nstokes, npix), outOftheta_val, dtype=float)
        for s in range(nstokes):
            beam_map[s, :len(theta)] = _get_beam_polar_value(beampolar, theta, phi, s)
        return BeamMap(beam_map)

    def plot(self, stokes="I", color_resol=20, figsize=6, cmap="jet", return_fields=False):
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

class BeamMap:
    """Represents a beam map.

    Attributes:
        nside (int): The resolution parameter of the map.
        map (numpy.ndarray): The beam map data.

    Methods:
        to_blm(lmax, mmax): Converts the beam map to spherical harmonic coefficients.

    """

    def __init__(self, map):
        self.nside = hp.npix2nside(map.shape[1])
        self.map = map

    def to_blm(self, lmax, mmax):
        """Converts the beam map to spherical harmonic coefficients.

        Args:
            lmax (int): The maximum multipole moment.
            mmax (int): The maximum value for a coefficient index.

        Returns:
            numpy.ndarray: The spherical harmonic coefficients.

        Raises:
            AssertionError: If lmax is greater than 3*nside-1 or if mmax is greater than lmax.

        """
        assert lmax <= 3*self.nside-1, "Error in BeamMap.to_blm: lmax < 3*nside-1"
        assert mmax <= lmax, "Error in BeamMap.to_blm: mmax < lmax"

        blm = hp.map2alm(self.map, lmax=lmax, mmax=mmax)
        return blm


def _get_beam_polar_value(beam:BeamPolar, theta, phi, s):
    """Calculate the value of the beam at a given theta, phi, and Stokes parameter.

    Args:
        beam (BeamPolar): The polar beam object.
        theta (float or array-like): The theta value(s) at which to evaluate the beam.
        phi (float or array-like): The phi value(s) at which to evaluate the beam.
        s (int): The Stokes parameter index.

    Returns:
        value (float or array-like): The value(s) of the beam at the given theta, phi, and Stokes parameter.

    """
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

def grasp2blm(
    filepath,
    nside,
    lmax,
    mmax,
    outOftheta_val=0.0,
    copol_axis='x',
    healpy_pol_convention=True
    ):
    """Convert a GRASP file to a spherical harmonic coefficients of beam map.

    Args:
        filepath (str): Path to the GRASP file.
        nside (int): Resolution parameter for the output beam map.
        lmax (int): Maximum l value for the spherical harmonic expansion.
        mmax (int): Maximum m value for the spherical harmonic expansion.
        outOftheta_val (float): Value to assign to pixels outside
            the valid theta range.
        copol_axis (str, optional): Axis of the co-polarization
            component. Defaults to 'x'.
        healpy_pol_convention (bool, optional): Whether to use the Healpix
            polarization convention. Defaults to True.

    Returns:
        blm (numpy.ndarray): Spherical harmonic coefficients of the beam map.

    Raises:
        ValueError: If the file format is unknown.

    """
    if filepath.endswith(".grd"):
        beam = BeamGrid(filepath)
    elif filepath.endswith(".cut"):
        beam = BeamCut(filepath)
    else:
        raise ValueError("Error in grasp2blm: unknown file format")
    beam_polar = beam.to_polar(copol_axis)
    beam_map = beam_polar.to_map(
        nside,
        healpy_pol_convention=healpy_pol_convention,
        outOftheta_val=outOftheta_val
        )
    blm = beam_map.to_blm(lmax, mmax)
    return blm
