import unittest
import os
from pathlib import Path
import numpy as np
import healpy as hp
import grasp2alm as g2a


class TestBeamCut(unittest.TestCase):
    """
    Testing the functionality of converting a gaussian beam
    to alm coefficients using the grasp2alm function.

    Methods:
        setUp():
            Initializes the test environment. This includes setting up parameters
            for the beam and writing the beam cut file.

        tearDown():
            Cleans up the test environment by removing the generated beam cut file.

        test_cut_beam2alm():
            Tests the grasp2alm function by comparing the generated alm coefficients
            from the beam cut file to the ideal alm coefficients of a Gaussian beam.
            Asserts that the coefficients match within a specified tolerance.
    """

    def setUp(self):
        """
        Sets up the test environment by initializing parameters.
        This includes defining the path to the beam cut file, setting the polarization flag,
        nside parameter, lmax parameter, and calculating the beam FWHM and sigma.
        It also generates the beam profile and writes it to a file.
        """
        self.path_cut: str = str(Path(__file__).parent / "beam_files" / "beam2alm.cut")
        self.path_grid: str = str(Path(__file__).parent / "beam_files" / "beam2alm.grd")

        self.nside: int = 1024
        self.pol: bool = True
        self.lmax: int = 2 * self.nside
        self.mmax: int = 2
        self.ellipticity: float = 1.0
        self.psi_ell_rad: float = 0.0
        self.psi_pol_rad: float = 0.0
        self.cross_polar_leakage: float = 0.0

        fwhm_deg: float = np.rad2deg(hp.nside2resol(self.nside)) * 50
        beam_sigma: float = np.deg2rad(fwhm_deg) / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        amplitude: float = 1 / (2 * np.pi * beam_sigma * beam_sigma)

        gauss = g2a.BeamGauss(fwhm_deg=fwhm_deg, amplitude=amplitude)
        gauss.write2cut(self.path_cut, -fwhm_deg * 2, 30001, 40)
        gauss.write2thetaphigrid(
            self.path_grid, 0.0, 0.0, 360.0, fwhm_deg * 2, 60, 30001
        )

        self.ideal_alm = gauss.get_alm(
            self.lmax,
            self.mmax,
            self.ellipticity,
            self.psi_ell_rad,
            self.psi_pol_rad,
            self.cross_polar_leakage,
        )

    def tearDown(self):
        """
        Cleans up the test environment by removing the generated beam cut file.
        """
        if os.path.exists(self.path_cut):
            os.remove(self.path_cut)
        if os.path.exists(self.path_grid):
            os.remove(self.path_grid)

    def test_cut_beam2alm(self):
        """
        Test function to verify if the alm coefficients generated from the beam cut file
        match the ideal alm coefficients of a Gaussian beam. The test asserts that the
        T, E, and B mode coefficients match within a specified tolerance of 1e-3.

        Asserts:
            True if the alm coefficients generated from the beam cut file
            are close to the ideal alm coefficients within a tolerance of 1e-3.
        """

        test_alm = g2a.grasp2alm(
            self.path_cut,
            self.nside,
            interp_method="pchip",
            copol_axis="y",
            lmax=self.lmax,
            mmax=self.mmax,
            pol=self.pol,
            use_pixel_weights=True,
        )

        # Get the indices for T, E and B mode coefficients
        index_T = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 0)
        index_E = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 2)
        index_B = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 2)

        # Assert that the alm coefficients match within a tolerance of 1e-3
        self.assertTrue(
            np.allclose(
                test_alm[0][
                    index_T
                ],  # T mode coefficients from the generated beam cut file
                self.ideal_alm[0][
                    index_T
                ],  # T mode coefficients from the ideal Gaussian beam
                atol=1e-3,  # Tolerance for the comparison
            )
        )
        self.assertTrue(
            np.allclose(
                test_alm[1][
                    index_E
                ],  # E mode coefficients from the generated beam cut file
                self.ideal_alm[1][
                    index_E
                ],  # E mode coefficients from the ideal Gaussian beam
                atol=1e-3,  # Tolerance for the comparison
            )
        )
        self.assertTrue(
            np.allclose(
                test_alm[2][
                    index_B
                ],  # B mode coefficients from the generated beam cut file
                self.ideal_alm[2][
                    index_B
                ],  # B mode coefficients from the ideal Gaussian beam
                atol=1e-3,  # Tolerance for the comparison
            )
        )

    def test_grid_beam2alm(self):
        """
        Test function to verify if the alm coefficients generated from the beam grid file
        match the ideal alm coefficients of a Gaussian beam. The test asserts that the
        T, E, and B mode coefficients match within a specified tolerance of 1e-3.

        Asserts:
            True if the alm coefficients generated from the beam grid file
            are close to the ideal alm coefficients within a tolerance of 1e-3.
        """

        test_alm = g2a.grasp2alm(
            self.path_grid,
            self.nside,
            interp_method="cubic",
            copol_axis="y",
            lmax=self.lmax,
            mmax=self.mmax,
            pol=self.pol,
            use_pixel_weights=True,
        )

        # Get the indices for T, E and B mode coefficients
        index_T = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 0)
        index_E = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 2)
        index_B = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 2)

        # Assert that the alm coefficients match within a tolerance of 1e-3
        self.assertTrue(
            np.allclose(
                test_alm[0][
                    index_T
                ],  # T mode coefficients from the generated beam grid file
                self.ideal_alm[0][
                    index_T
                ],  # T mode coefficients from the ideal Gaussian beam
                atol=1e-3,  # Tolerance for the comparison
            )
        )
        self.assertTrue(
            np.allclose(
                test_alm[1][
                    index_E
                ],  # E mode coefficients from the generated beam grid file
                self.ideal_alm[1][
                    index_E
                ],  # E mode coefficients from the ideal Gaussian beam
                atol=1e-3,  # Tolerance for the comparison
            )
        )
        self.assertTrue(
            np.allclose(
                test_alm[2][
                    index_B
                ],  # B mode coefficients from the generated beam grid file
                self.ideal_alm[2][
                    index_B
                ],  # B mode coefficients from the ideal Gaussian beam
                atol=1e-3,  # Tolerance for the comparison
            )
        )


if __name__ == "__main__":
    unittest.main()
