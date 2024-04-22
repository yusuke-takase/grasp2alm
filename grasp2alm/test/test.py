import unittest
import grasp2alm as gp
import numpy as np

class TestGrasp2Alm(unittest.TestCase):
    def test_grasp2alm_with_gridfile(self):
        file = "./beam_files/test.grd"
        beam = gp.BeamGrid(file)
        polar = beam.to_polar()
        nside = 128
        beammap = polar.to_map(nside)
        blm1 = beammap.to_alm()
        blm2 = gp.grasp2alm(file, nside)
        self.assertTrue(np.array_equal(blm1, blm2), "blm1 and blm2 are not equal")

    def test_grasp2alm_with_cutfile(self):
        file = "./beam_files/test.cut"
        beam = gp.BeamCut(file)
        polar = beam.to_polar()
        nside = 128
        beammap = polar.to_map(nside)
        blm1 = beammap.to_alm()
        blm2 = gp.grasp2alm(file, nside)
        self.assertTrue(np.array_equal(blm1, blm2), "blm1 and blm2 are not equal")


if __name__ == '__main__':
    unittest.main()
