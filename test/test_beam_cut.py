import unittest
import numpy as np
from pathlib import Path
from grasp2alm import BeamCut

class TestBeamCut(unittest.TestCase):
    def setUp(self):
        path_to_test_cut = str(Path(__file__).parent / "beam_files" / "unit_test.cut")
        txt = '''Test header
            -180 90 3 0 3 1 2
            1 2 3 4
            -1.0 -2.0 -3.0 -4.0
            1 1 1 1'''
        text_file = open(path_to_test_cut, "w")
        text_file.write(txt)
        text_file.close()
        self.test_cut = BeamCut(path_to_test_cut)
    def test_cut_reading(self):
        assert self.test_cut.header == "Test header"
        assert self.test_cut.vini == -180.0
        assert self.test_cut.vinc == 90.0
        assert self.test_cut.vnum == 3
        assert self.test_cut.c == np.array([0.0])
        assert self.test_cut.icomp == 3
        assert self.test_cut.icut == 1
        assert self.test_cut.ncomp == 2
        assert self.test_cut.ncut == 1
        self.assertTrue(
            np.array_equal(
                self.test_cut.amp,
                np.array(
                    [[[1.+2.j],[-1.-2.j],[1.+1.j]],
                    [[3.+4.j],[-3.-4.j],[1.+1.j]]]
                )
            )
        )

if __name__ == '__main__':
    unittest.main()