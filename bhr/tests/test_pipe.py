import unittest

from bhr.pipe import Pipe


class TestPipe(unittest.TestCase):

    def test_init_ip_pipe(self):
        inputs = {
            "fluid-type": "WATER",
            "nominal-pipe-diameter-inch": 1.0,
            "length": 100
        }

        pipe = Pipe(inputs)

        self.assertEqual(pipe.length, 100)
        self.assertAlmostEqual(pipe.outer_diameter, 0.03340, delta=0.0001)
        self.assertAlmostEqual(pipe.inner_diameter, 0.02732, delta=0.0001)
        self.assertAlmostEqual(pipe.area_s_outer, 10.49, delta=0.01)
        self.assertAlmostEqual(pipe.area_s_inner, 8.58, delta=0.01)

    def test_calc_resist(self):
        inputs = {
            "fluid-type": "WATER",
            "nominal-pipe-diameter-inch": 1.0,
            "length": 100
        }

        pipe = Pipe(inputs)

        self.assertAlmostEqual(pipe.calc_resist(0.5, 20), 0.0906, delta=0.0001)
