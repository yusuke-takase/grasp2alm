import unittest
import numpy as np
from pathlib import Path
from grasp2alm import BeamCut

class TestBeamCut(unittest.TestCase):
    def setUp(self):
        self.path_to_test_cut = str(Path(__file__).parent / "beam_files" / "unit_test.cut")
    def write_to_test_cut(self, txt:str):
        text_file = open(self.path_to_test_cut, "w")
        text_file.write(txt)
        text_file.close()
    def test_input_extension_exeption(self):
        with self.assertRaises(ValueError):
            BeamCut("test.not_cut")
    def test_ncomp_exeption(self):
        txt_with_error = '''Test header
            -180 90 3 0 3 1 100
            1 2 3 4
            -1.0 -2.0 -3.0 -4.0
            1 1 1 1
        '''
        self.write_to_test_cut(txt_with_error)
        with self.assertRaises(ValueError):
            BeamCut(self.path_to_test_cut)
    def test_vnum_exeption(self):
        txt_with_error = '''Test header
            -180 90 10 0 3 1 2
            1 2 3 4
            -1.0 -2.0 -3.0 -4.0
            1 1 1 1
        '''
        self.write_to_test_cut(txt_with_error)
        with self.assertRaises(ValueError):
            BeamCut(self.path_to_test_cut)
    def test_cut_reading(self):
        txt = '''Test header
            -180 90 3 0 3 1 2
            1 2 3 4
            -1.0 -2.0 -3.0 -4.0
            1 1 1 1
        '''
        self.write_to_test_cut(txt)
        test_cut = BeamCut(self.path_to_test_cut)
        assert test_cut.header == "Test header"
        assert test_cut.vini == -180.0
        assert test_cut.vinc == 90.0
        assert test_cut.vnum == 3
        assert test_cut.c == np.array([0.0])
        assert test_cut.icomp == 3
        assert test_cut.icut == 1
        assert test_cut.ncomp == 2
        assert test_cut.ncut == 1
        expected_amp = np.array(
                    [[[1.+2.j],[-1.-2.j],[1.+1.j]],
                    [[3.+4.j],[-3.-4.j],[1.+1.j]]]
                )
        self.assertTrue(np.array_equal(test_cut.amp,expected_amp))

if __name__ == '__main__':
    unittest.main()