from unittest import TestCase

from bhr.u_tube import UTube


class TestUTube(TestCase):
    def test_init(self):
        inputs = {
            "pipe_outer_diameter": 0.034,
            "pipe_dimension_ratio": 11,
            "length": 100,
            "shank_space": 0.03,
            "pipe_conductivity": 0.4,
            "fluid_type": "WATER",
        }

        u_tube = UTube(**inputs)

        self.assertEqual(u_tube.length, 100)
        self.assertEqual(u_tube.pipe_length, 200)
        self.assertEqual(u_tube.shank_space, 0.03)
