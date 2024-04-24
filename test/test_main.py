# -*- encoding: utf-8 -*-

import unittest
from grasp2alm import BeamGrid, BeamCut, grasp2alm
import numpy as np
from pathlib import Path

class TestGrasp2Alm(unittest.TestCase):
    def test_grasp2alm_with_gridfile(self):
        file = Path(__file__).parent / "beam_files" / "test.grd"
        file = str(file)
        beam = BeamGrid(file)
        polar = beam.to_polar()
        nside = 128
        beammap = polar.to_map(nside)
        blm1 = beammap.to_alm()
        blm2 = grasp2alm(file, nside)
        self.assertTrue(np.array_equal(blm1, blm2), "blm1 and blm2 are not equal")

    def test_grasp2alm_with_cutfile(self):
        file = Path(__file__).parent / "beam_files" / "test.cut"
        file = str(file)
        beam = BeamCut(file)
        polar = beam.to_polar()
        nside = 128
        beammap = polar.to_map(nside)
        blm1 = beammap.to_alm()
        blm2 = grasp2alm(file, nside)
        self.assertTrue(np.array_equal(blm1, blm2), "blm1 and blm2 are not equal")


if __name__ == '__main__':
    unittest.main()
