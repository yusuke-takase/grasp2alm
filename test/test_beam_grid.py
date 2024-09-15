import unittest
import os
from pathlib import Path
import numpy as np
from grasp2alm import BeamGrid

class TestBeamGrid(unittest.TestCase):
    """
    Unit tests for the BeamGrid class to ensure proper handling of beam grid files and
    correct exception raising for invalid inputs.
    
    Methods:
        setUp():
            Sets up the test environment by defining the path to the test grid file.
    
        tearDown():
            Cleans up the test environment by removing the test grid file after each test.
    
        write_to_test_grid(txt: str):
            Writes the provided text to the test grid file.
    
        test_input_extension_exception():
            Tests that the BeamGrid class raises a ValueError 
            for files with an incorrect extension.
    
        test_input_grid_format():
            Tests that the BeamGrid class raises a ValueError 
            for an incorrect number for grid format (ktype) in the grid file.
    
        test_input_beams_number():
            Tests that the BeamGrid class raises a Warning 
            for an incorrect number of beams (nset) in the grid file.
        
        test_input_beam_solid_angle():
            Tests that the BeamGrid class raises a Warning 
            for an incorrect value of beam solid angle in the grid file.

    
        test_nan_exception():
            Tests that the BeamGrid class raises a ValueError 
            when the grid file contains NaN values.
    
        test_grid_reading():
            Tests that the BeamGrid class correctly reads a 
            properly formatted grid file and verifies the contents.
"""

    def setUp(self):
        """
        Sets up the test environment by defining the path to the test grid file.
        """
        self.path_to_test_grid = str(Path(__file__).parent / "beam_files" / "unit_test.grd")

    def tearDown(self):
        """
        Cleans up the test environment by removing 
        the test grid file after each test method is executed.
        """
        if os.path.exists(self.path_to_test_grid):
            os.remove(self.path_to_test_grid)

    def write_to_test_grid(self, txt:str):
        """
        Writes the provided text to the test grid file.

        Args:
            txt (str): The text to write to the test cut file.
        """
        text_file = open(self.path_to_test_grid, "w", encoding='utf-8')
        text_file.write(txt)
        text_file.close()

    def test_input_extension_exception(self):
        '''
        Tests that the BeamGrid class raises a ValueError for files with an incorrect extension.

        Asserts:
            Raises ValueError when the file extension is not .grid.
        '''
        with self.assertRaises(ValueError):
            BeamGrid("test.not_grid")

    def test_input_grid_format(self):
        '''
        Tests that the BeamGrid class raises a ValueError for an incorrect number for grid format in the grid file.

        Asserts:
            Raises ValueError when the grid format (ktype) is not 1.
        '''
        txt_with_error = "Test header\n++++\n2"
        self.write_to_test_grid(txt_with_error)
        with self.assertRaises(ValueError):
            BeamGrid(self.path_to_test_grid)

    def test_input_beams_number(self):
        '''
        Tests that the BeamGrid class raises a Warning for an incorrect number of beams in the grid file.
        
        Asserts:
            Raises a Warning when the number of beams (nset) is not 1.
        '''
        txt_with_error = \
            "Test header\n" + \
            "++++\n" + \
            "1\n" + \
            "2 3 2 7\n" + \
            "0 0\n\n" + \
            "0.0 0.0 360.0 90.0\n" + \
            "2 2 0\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1"
        self.write_to_test_grid(txt_with_error)
        with self.assertWarns(Warning):
            BeamGrid(self.path_to_test_grid)

    def test_input_beam_solid_angle(self):
        '''
        Tests that the BeamGrid class raises a Warning for an incorrect value of beam solid angle in the grid file.

        Asserts:
            Raises a Warning if the beam solid angle is different from 2pi or 4pi
        '''
        txt_with_error = \
            "Test header\n" + \
            "++++\n" + \
            "1\n" + \
            "1 3 2 7\n" + \
            "0 0\n" + \
            "0.0 0.0 340.0 80.0\n" + \
            "2 2 0\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1"        
        self.write_to_test_grid(txt_with_error)
        with self.assertWarns(Warning):
            BeamGrid(self.path_to_test_grid)

    def test_nan_exception(self):
        '''
        Tests that the BeamGrid class raises a ValueError 
        when the grid file contains NaN values.

        Asserts:
            Raises a ValueError when the gridfile contains NaN values

        '''
        txt_with_error = \
            "Test header\n" + \
            "++++\n" + \
            "1\n" + \
            "1 3 2 7\n" + \
            "0 0\n" + \
            '0.0 0.0 360.0 90.0\n' + \
            "2 2 0\n" + \
            "1 1 1 1\n" + \
            "1 Nan 1 1\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1"
        self.write_to_test_grid(txt_with_error)
        with self.assertRaises(ValueError):
            BeamGrid(self.path_to_test_grid)

    def test_grid_reading(self):
        '''
        Tests that the BeamGrid class correctly reads a 
        properly formatted grid file and verifies the contents.

        Asserts:
            Asserts that the header, ktype, nset, klimit,
            icomp, ncomp, igrid, ix, iy, xs, ys, xe, 
            ye, nx, ny, freq, frequnit
            are correctly read from the grid file
        '''
        txt = \
        "Test header\n" + \
        "VERSION: TICRA-EM-FIELD-V0.1\n" + \
        "FREQUENCY_NAME: freq\n" + \
        "FREQUENCIES [GHz]:\n" + \
        "119.0\n" + \
        "++++\n" + \
        "1\n" + \
        "1 3 2 7\n" + \
        "0 0\n" + \
        '0.0 0.0 360.0 90.0\n' + \
        "3 3 0\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1"

        self.write_to_test_grid(txt)
        test_grid = BeamGrid(self.path_to_test_grid)

        assert test_grid.header == "Test header\n" + \
                                  "VERSION: TICRA-EM-FIELD-V0.1\n" + \
                                  "FREQUENCY_NAME: freq\n"
        assert test_grid.ktype == 1
        assert test_grid.nset == 1
        assert test_grid.klimit == 0
        assert test_grid.icomp == 3
        assert test_grid.ncomp == 2
        assert test_grid.igrid == 7
        assert test_grid.ix == 0
        assert test_grid.iy == 0
        assert test_grid.xs == 0.0
        assert test_grid.ys == 0.0
        assert test_grid.xe == 360.0
        assert test_grid.ye == 90.0
        assert test_grid.nx == 3
        assert test_grid.ny == 3
        assert test_grid.freq == 119.0
        assert test_grid.frequnit == "GHz"

        expected_amp = np.array([
            [[1.+1.j,1.+1.j,1.+1.j],
            [1.+1.j,1.+1.j,1.+1.j],
            [1.+1.j,1.+1.j,1.+1.j]],
  
            [[1.+1.j,1.+1.j,1.+1.j],
            [1.+1.j,1.+1.j,1.+1.j],
            [1.+1.j,1.+1.j,1.+1.j]],
        ])
        self.assertTrue(np.array_equal(test_grid.amp,expected_amp))

if __name__ == "__main__":
    unittest.main()
