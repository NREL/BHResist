import unittest

from bhr.pipe import Pipe


class TestPipe(unittest.TestCase):

    def test_init_pipe(self):
        inputs = {
            "pipe_outer_diameter": 0.0334,
            "pipe_dimension_ratio": 11,
            "pipe_length": 100,
            "pipe_conductivity": 0.4,
            "fluid_type": "WATER"
        }

        pipe = Pipe(**inputs)

        self.assertEqual(pipe.pipe_length, 100)
        self.assertAlmostEqual(pipe.pipe_outer_diameter, 0.03340, delta=0.0001)
        self.assertAlmostEqual(pipe.pipe_inner_diameter, 0.02732, delta=0.0001)
        self.assertAlmostEqual(pipe.area_s_outer, 10.49, delta=0.01)
        self.assertAlmostEqual(pipe.area_s_inner, 8.58, delta=0.01)

    def test_calc_resist(self):
        inputs = {
            "pipe_outer_diameter": 0.0334,
            "pipe_dimension_ratio": 11,
            "pipe_length": 100,
            "pipe_conductivity": 0.4,
            "fluid_type": "WATER"
        }

        pipe = Pipe(**inputs)

        self.assertAlmostEqual(pipe.calc_pipe_resist(0.5, 20), 0.0906, delta=0.0001)
