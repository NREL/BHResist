from unittest import TestCase

from bhr.coaxial_borehole import Coaxial


class TestCoaxialBorehole(TestCase):
    def setUp(self):
        self.inputs = {
                 "borehole_diameter": float,
                 "outer_pipe_outer_diameter": float,
                 "outer_pipe_dimension_ratio": float,
                 "outer_pipe_conductivity": float,
                 "inner_pipe_outer_diameter": float,
                 "inner_pipe_dimension_ratio": float,
                 "inner_pipe_conductivity": float,
                 "length": 100,
                 "grout_conductivity": float,
                 "soil_conductivity": float,
                 "fluid_type": str,
                 "fluid_concentration": float,
        }

    def test_init(self):
        bh = Coaxial(**self.inputs)
        self.assertEqual(bh.length, 100)