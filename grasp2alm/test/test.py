import unittest
import grasp2alm as gp
import numpy as np

#gridfile = "../../../notebook/pix0000_119_pp_f2p2_v4_mft_tp.grd"
#cutfile = "../../../notebook/MFT_88.5GHz_000.0_166.7_xpol_v2.cut"

class TestGrasp2Alm(unittest.TestCase):
    # TODO : add dummy file in the ./beam_files
    def test_grasp2alm_with_gridfile(self):
        file = "../../../notebook/pix0000_119_pp_f2p2_v4_mft_tp.grd"
        beam = gp.BeamGrid(file)
        polar = beam.to_polar()
        nside = 128
        beammap = polar.to_map(nside)
        blm1 = beammap.to_alm()
        blm2 = gp.grasp2alm(file, nside)
        self.assertTrue(np.array_equal(blm1, blm2), "blm1 and blm2 are not equal")

    def test_grasp2alm_with_cutfile(self):
        file = "../../../notebook/MFT_88.5GHz_000.0_166.7_xpol_v2.cut"
        beam = gp.BeamCut(file)
        polar = beam.to_polar()
        nside = 128
        beammap = polar.to_map(nside)
        blm1 = beammap.to_alm()
        blm2 = gp.grasp2alm(file, nside)
        self.assertTrue(np.array_equal(blm1, blm2), "blm1 and blm2 are not equal")


if __name__ == '__main__':
    unittest.main()
