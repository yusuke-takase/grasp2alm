import unittest
import os
from pathlib import Path
import numpy as np
from grasp2alm import BeamCut

class TestBeamCut(unittest.TestCase):
    """
    Unit tests for the BeamCut class to ensure proper handling of beam cut files and
    correct exception raising for invalid inputs.

    Methods:
        setUp():
            Sets up the test environment by defining the path to the test cut file.
        
        tearDown():
            Cleans up the test environment by removing the test cut file after each test.
        
        write_to_test_cut(txt: str):
            Writes the provided text to the test cut file.
        
        test_input_extension_exception():
            Tests that the BeamCut class raises a ValueError 
            for files with an incorrect extension.
        
        test_ncomp_exception():
            Tests that the BeamCut class raises a ValueError 
            for an incorrect number of components in the cut file.
        
        test_vnum_exception():
            Tests that the BeamCut class raises a ValueError 
            for an incorrect number of data points (vnum) in the cut file.
        
        test_nan_exception():
            Tests that the BeamCut class raises a ValueError 
            when the cut file contains NaN values.
        
        test_cut_reading():
            Tests that the BeamCut class correctly reads a 
            properly formatted cut file and verifies the contents.
    """

    def setUp(self):
        """
        Sets up the test environment by defining the path to the test cut file.
        """
        self.path_to_test_cut = str(Path(__file__).parent / "beam_files" / "unit_test.cut")

    def tearDown(self):
        """
        Cleans up the test environment by removing 
        the test cut file after each test method is executed.
        """
        if os.path.exists(self.path_to_test_cut):
            os. remove(self.path_to_test_cut)

    def write_to_test_cut(self, txt:str):
        """
        Writes the provided text to the test cut file.

        Args:
            txt (str): The text to write to the test cut file.
        """
        with open(self.path_to_test_cut, "w", encoding='utf-8') as text_file:
            text_file.write(txt)

    def test_input_extension_exeption(self):
        """
        Tests that the BeamCut class raises a ValueError for files with an incorrect extension.

        Asserts:
            Raises ValueError when the file extension is not .cut.
        """
        with self.assertRaises(ValueError):
            BeamCut("test.not_cut")

    def test_ncomp_exeption(self):
        """
        Tests that the BeamCut class raises a ValueError 
        for an incorrect number of components in the cut file.

        Asserts:
            Raises ValueError when the number of components (ncomp) is incorrect.
        """
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
        """
        Tests that the BeamCut class raises a ValueError 
        for an incorrect number of data points (vnum) in the cut file.

        Asserts:
            Raises ValueError when the number of data points (vnum) is incorrect.
        """
        txt_with_error = '''Test header
            -180 90 10 0 3 1 2
            1 2 3 4
            -1.0 -2.0 -3.0 -4.0
            1 1 1 1
        '''
        self.write_to_test_cut(txt_with_error)
        with self.assertRaises(ValueError):
            BeamCut(self.path_to_test_cut)

    def test_nan_exeption(self):
        """
        Tests that the BeamCut class raises a ValueError 
        when the cut file contains NaN values.

        Asserts:
            Raises ValueError when the cut file contains NaN values.
        """
        txt_with_error = '''Test header
            -180 90 3 0 3 1 2
            1 2 3 4
            -1.0 nan -3.0 -4.0
            1 1 1 1
        '''
        self.write_to_test_cut(txt_with_error)
        with self.assertRaises(ValueError):
            BeamCut(self.path_to_test_cut)

    def test_cut_reading(self):
        """
        Tests that the BeamCut class correctly reads a properly formatted cut file
        and verifies the contents.

        Asserts:
            Asserts that the header, vini, vinc, vnum, c, 
            icomp, icut, ncomp, ncut, and amp
            are correctly read from the cut file.
        """
        txt = '''Test header
            -180 90 3 0 3 1 2
            1 2 3 4
            -1.0 -2.0 -3.0 -4.0
            1 1 1 1
        '''
        self.write_to_test_cut(txt)
        test_cut = BeamCut(self.path_to_test_cut)
        expected_amp = np.array(
                    [[[1.+2.j],[-1.-2.j],[1.+1.j]],
                    [[3.+4.j],[-3.-4.j],[1.+1.j]]]
                )
        self.assertTrue(test_cut.header == "Test header")
        self.assertTrue(test_cut.vini == -180.0)
        self.assertTrue(test_cut.vinc == 90.0)
        self.assertTrue(test_cut.vnum == 3)
        self.assertTrue(test_cut.c == np.array([0.0]))
        self.assertTrue(test_cut.icomp == 3)
        self.assertTrue(test_cut.icut == 1)
        self.assertTrue(test_cut.ncomp == 2)
        self.assertTrue(test_cut.ncut == 1)
        self.assertTrue(np.array_equal(test_cut.amp,expected_amp))

if __name__ == '__main__':
    unittest.main()
