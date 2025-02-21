import numpy as np
import healpy as hp
from healpy import Alm
from typing import Optional, Union
from scipy.special import iv as bessel_i


def alm_index(lmax: int, ell: int, m: int) -> int:
    """Return the index of an a_ℓm coefficient in an array"""

    return m * (2 * lmax + 1 - m) // 2 + ell


class BeamGauss:
    """Class to generate a gaussian beam."""

    def __init__(self, fwhm_deg: float, amplitude: float = 1.0):
        """Initializes a BeamGauss object.

        Args:
            amplitude (`float`): Amplitude of the beam.
            fwhm (`float`): Full width at half maximum (FWHM) of the beam.

        Raises:
            ValueError: If amplitude or fwhm is not a number.
        """

        assert (
            isinstance(amplitude, (int, float)) is True
        ), "Amplitude must be a number."
        assert isinstance(fwhm_deg, (int, float)) is True, "FWHM must be a number."

        self.amplitude = amplitude
        self.fwhm_deg = fwhm_deg

    def gaussian_beam(self, theta: float):
        """Calculates the value of a Gaussian beam at a given angle.

        Args:
            theta (`float`): The angle at which to calculate the beam value.

        Returns:
            `float`: The value of the Gaussian beam at the given angle.
        """
        fwhm_rad = np.deg2rad(self.fwhm_deg)
        sigmasq = fwhm_rad * fwhm_rad / (8 * np.log(2.0))
        return self.amplitude * np.exp(-0.5 * theta**2 / sigmasq)

    def write2cut(self, path: str, vini: int, vnum: int, ncut: int):
        """Writes the formatted beam data to the specified path with the provided headers.
        Each cut contains vnum lines of data.

        Args:
            path (`str`): The path to write the beam data to.
            vnum (`int`): The number of lines in each cut.
            ncut (`int`): The number of cuts.
            co (array): The beam data.

        Raises:
            ValueError: If the path is not a string.
        """
        vinc: float = abs(vini) * 2 / vnum
        header_1: str = "Field data in cuts"

        theta = np.linspace(vini, -vini, vnum)
        theta = np.deg2rad(theta)
        phi = np.linspace(0, 180 - 180 / ncut, ncut)
        beam = self.gaussian_beam(theta)

        assert isinstance(path, str) is True, "Path must be a string."
        with open(path, "w", encoding="utf-8") as file:
            for j in range(ncut):
                file.write(header_1)
                file.write("\n")
                file.write(f"{vini} {vinc} {vnum} {phi[j]} 3 1 2")
                file.write("\n")
                for i in range(vnum):
                    co_i = np.emath.sqrt(beam[i])
                    cx_i = 0.0
                    file.write(
                        f"{np.real(co_i)} {np.imag(co_i)} {np.real(cx_i)} {np.imag(cx_i)}\n"
                    )

    def write2thetaphigrid(
        self, path: str, xs: float, ys: float, xe: float, ye: float, nx: int, ny: int
    ):
        """Writes the formatted beam data to the specified path with the provided headers.
        Each cut contains vnum lines of data.

        Args:
            path (`str`): The path to write the beam data to.
            xs (`float`): The starting x-coordinate of the grid.
            ys (`float`): The starting y-coordinate of the grid.
            xe (`float`): The ending x-coordinate of the grid.
            ye (`float`): The ending y-coordinate of the grid.
            nx (`int`): The number of grid points in the x-direction.
            ny (`int`): The number of grid points in the y-direction.
        """
        header: str = (
            "VERSION: TICRA-EM-FIELD-V0.1\n"
            + "Field data in grid\n"
            + "SOURCE_FIELD_NAME: baffle_aperture.po\n"
            + "FREQUENCY_NAME: freq\n"
            + "FREQUENCIES [GHz]:\n"
            + "30\n"
            + "++++\n"
            + "1\n"
            + "1 3 2 7\n"
            + "0 0\n"
            + f"{xs} {ys} {xe} {ye}\n"
            + f"{nx} {ny} 0"
        )
        phi = np.linspace(xs, xe, nx, endpoint=False)
        theta = np.linspace(ys, ye, ny, endpoint=False)
        grid = np.deg2rad(np.meshgrid(phi, theta))
        beam = self.gaussian_beam(grid[1])

        assert isinstance(path, str) is True, "Path must be a string."
        with open(path, "w", encoding="utf-8") as file:
            file.write(header)
            file.write("\n")
            co = np.emath.sqrt(beam.flatten())
            for i in co:
                file.write(f"{np.real(i)} {np.imag(i)} 0.0 0.0\n")
        return grid

    def get_alm(
        self,
        lmax: int,
        mmax: int,
        ellipticity: Optional[float] = 1.0,
        psi_ell_rad: Optional[float] = 0.0,
        psi_pol_rad: Union[float, None] = 0.0,
        cross_polar_leakage: Optional[float] = 0.0,
    ) -> np.ndarray:
        """Return an array of spherical harmonics :math:`a_{\ell m}` that represents a Gaussian beam

        Args:
            lmax (`int`): The maximum value for :math:`\ell`
            mmax (`int`): The maximum range for :math:`m`; usually this is equal to ``lmax``
            ellipticity (`float`): The ellipticity of the beam, defined as major axis/minor
            axis. Default is 1 (circular beam).
            psi_ell_rad (`float`): The inclination of the major axis of the ellipse with
            respect to the x-axis. This is not relevant for cirular beams. Default is 0.
            psi_pol_rad (`float`): The polarization of the beam with respect to the x-axis,
            if `None` only I beam will be returned. Default is 0 rad.
            cross_polar_leakage (`float`): The cross-polar leakage (pure number). Default is 0.

        Return:
            :math:`a_{\ell m}` values (`numpy.ndarray`)
        """
        is_elliptical = False if ellipticity == 1.0 else True

        is_polarized = psi_pol_rad is not None

        if is_polarized:
            num_stokes = 3
        else:
            num_stokes = 1

        nval = Alm.getsize(lmax, mmax)
        alm = np.zeros((num_stokes, nval), dtype=np.complex128)
        fwhm_rad = np.deg2rad(self.fwhm_deg)

        if not is_elliptical:
            # Circular beam
            sigma_squared = fwhm_rad**2 / (8 * np.log(2))
            for ell in range(lmax + 1):
                alm[0, alm_index(lmax, ell, 0)] = np.sqrt(
                    (2 * ell + 1) / (4 * np.pi)
                ) * np.exp(-0.5 * sigma_squared * ell * (ell + 1))

            if is_polarized:
                f1 = np.cos(2 * psi_pol_rad) - np.sin(2 * psi_pol_rad) * 1.0j
                for ell in range(2, lmax + 1):
                    value = (
                        np.sqrt((2 * ell + 1) / (32 * np.pi))
                        * np.exp(-0.5 * sigma_squared * ell * (ell + 1))
                        * f1
                    )
                    alm[1, alm_index(lmax, ell, 2)] = value
                    alm[2, alm_index(lmax, ell, 2)] = value * 1.0j
        else:
            # Elliptical beam
            e_squared = 1.0 - 1.0 / ellipticity**2
            sigma_x_squared = fwhm_rad**2 * ellipticity / (np.log(2) * 8)

            # I component
            for ell in range(lmax + 1):
                tmp = ell * (ell + 1) * sigma_x_squared
                for m in range(0, min(ell, mmax) + 1, 2):
                    alm[0, alm_index(lmax, ell, m)] = (
                        np.sqrt((2 * ell + 1) / (4 * np.pi))
                        * np.exp(-0.5 * tmp * (1.0 - e_squared / 2))
                        * bessel_i(m // 2, 0.25 * tmp * e_squared)
                    )

            if is_polarized:
                # Do G and C components

                rho = psi_pol_rad - psi_ell_rad
                f1 = np.cos(2 * rho) - np.sin(2 * rho) * 1j
                f2 = np.cos(2 * rho) + np.sin(2 * rho) * 1j

                for ell in range(2, lmax + 1):
                    tmp = ell * (ell + 1) * sigma_x_squared
                    tmp2 = 0.25 * tmp * e_squared

                    # m = 0
                    value = (
                        np.sqrt((2 * ell + 1) / (8 * np.pi))
                        * np.exp(-tmp * (1.0 - e_squared / 2) / 2)
                        * bessel_i(1, tmp2)
                    )

                    alm[1, alm_index(lmax, ell, 0)] = value * np.cos(2 * rho)
                    alm[2, alm_index(lmax, ell, 0)] = value * np.sin(2 * rho)

                    # m = 2, 4, …

                    for m in range(2, min(ell, mmax) + 1, 2):
                        value = np.sqrt((2 * ell + 1) / (8 * (4 * np.pi))) * np.exp(
                            -0.5 * tmp * (1.0 - 0.5 * e_squared)
                        )
                        b1 = f1 * bessel_i((m - 2) // 2, tmp2)
                        b2 = f2 * bessel_i((m + 2) // 2, tmp2)
                        alm[1, alm_index(lmax, ell, m)] = value * (b1 + b2)
                        alm[2, alm_index(lmax, ell, m)] = value * (b1 - b2) * 1j

            # Rotate the multipoles through angle psi_ell about the z-axis, so
            # the beam is in the right orientation (only need this for even m).

            for m in range(0, mmax + 1, 2):
                f1 = np.cos(m * psi_ell_rad) - np.sin(m * psi_ell_rad) * 1j
                for n in range(num_stokes):
                    for ell in range(m, lmax + 1):
                        alm[n, alm_index(lmax, ell, m)] *= f1

        # Adjust multipoles for cross-polar leakage
        alm[0, :] *= 1.0 + cross_polar_leakage

        if num_stokes > 1:
            for n in (1, 2):
                alm[n, :] *= 1.0 - cross_polar_leakage

        # Adjust the normalization
        if num_stokes > 1:
            for n in (1, 2):
                alm[n, :] *= -np.sqrt(2.0)
        self.alm = alm
        return alm

    def get_profile(self, nside: int):
        """Provide the beam map profile at a given nside.

        Args:
            nside (`int`): The nside of the map.

        Returns:
            `numpy.ndarray`: The beam profile map.
        """
        assert hasattr(self, "alm"), "Please run BeamGauss.get_alm() first"
        return hp.alm2map(self.alm, nside)
