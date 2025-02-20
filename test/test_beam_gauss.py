import unittest
from pathlib import Path
import numpy as np
from grasp2alm import BeamGauss


class TestBeamCut(unittest.TestCase):
    """
    Testing the functionality of converting a gaussian beam
    to alm coefficients using the grasp2alm function.

    Methods:
        setUp():
            Initializes the test environment. This includes setting up parameters
            for the beam and writing the beam cut file.

        test_get_alm():
            Tests the get_alm method of the BeamGauss class.
    """

    def setUp(self):
        """
        Sets up the parameters for the Gaussian beam.
        """
        self.path_gauss_alm: str = str(
            Path(__file__).parent / "reference" / "gauss_alm.npy"
        )
        self.lmax = 500
        self.mmax = self.lmax
        self.fwhm_deg = 1.0

        self.ellipticity = 0.5
        self.psi_ell_rad = np.pi / 6
        self.psi_pol_rad = np.pi / 4
        self.cross_polar_leakage = 0.0

    def test_get_alm(self, save: bool = False):
        """
        Tests the get_alm method of the BeamGauss class.

        Args:
            save (bool): If True, the generated alm coefficients are saved to a file.
        """
        Gauss = BeamGauss(self.fwhm_deg)
        alm = Gauss.get_alm(
            self.lmax,
            self.mmax,
            self.ellipticity,
            self.psi_ell_rad,
            self.psi_pol_rad,
            self.cross_polar_leakage,
        )
        if save:
            np.save(self.path_gauss_alm, alm)
        else:
            alm_ref = np.load(self.path_gauss_alm)
            self.assertTrue(
                np.allclose(alm, alm_ref),
                "Generated alm coefficients do not match reference values.",
            )
