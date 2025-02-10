import unittest

from bhr.borehole import Borehole


class TestBorehole(unittest.TestCase):

    def test_init_single_u_from_dict(self):
        inputs = {
            "fluid_type": "PROPYLENEGLYCOL",
            "fluid_concentration": 0.2,
            "boundary_condition": "uniform_heat_flux",
            "borehole_type": "single_u_tube",
            "single_u_tube": {
                "pipe_outer_diameter": 0.042,
                "pipe_dimension_ratio": 11,
                "pipe_conductivity": 0.4,
                "shank_space": 0.02,
            },
   
            "grout_conductivity": 1.2,
            "soil_conductivity": 2.5,
            "length": 100,
            "borehole_diameter": 0.14
        }

        bh = Borehole()
        bh.init_from_dict(inputs)

        # only pass flow rate, so pipe resistance should be computed in the process of this call

        self.assertAlmostEqual(bh.calc_bh_resist(temperature=20, flow_rate=0.5), 0.20425, delta=0.0001)

    def test_init_double_u_from_dict(self):
        inputs = {
            "fluid_type": "WATER",
            "fluid_concentration": 0,
            "boundary_condition": "uniform_heat_flux",
            "borehole_type": "double_u_tube",
            "double_u_tube": {
                "pipe_outer_diameter": 0.032,
                "pipe_dimension_ratio": 18.9,
                "pipe_conductivity": 0.389,
                "shank_space": 0.032,
                "pipe_inlet_arrangement": "ADJACENT"  # or DIAGONAL
            },
            "grout_conductivity": 1.5,
            "soil_conductivity": 3,
            "length": 200,
            "borehole_diameter": 0.115
        }

        bh = Borehole()
        bh.init_from_dict(inputs)

        # only pass flow rate, so pipe resistance should be computed in the process of this call
        self.assertAlmostEqual(bh.calc_bh_resist(temperature=20, flow_rate=0.2077), 0.1090, delta=0.0001)

    def test_init_coaxial_from_dict(self):
        inputs = {
            "fluid_type": "WATER",
            "fluid_concentration": 0,
            "boundary_condition": "uniform_heat_flux",
            "borehole_type": "coaxial",
            "coaxial": {
                "outer_pipe_outer_diameter": 0.064,
                "outer_pipe_dimension_ratio": 11,
                "outer_pipe_conductivity": 0.389,
                "inner_pipe_outer_diameter": 0.032,
                "inner_pipe_dimension_ratio": 11,
                "inner_pipe_conductivity": 0.389
            },
            "grout_conductivity": 1.5,
            "soil_conductivity": 3,
            "length": 200,
            "borehole_diameter": 0.115
        }

        bh = Borehole()
        bh.init_from_dict(inputs)

        # only pass flow rate, so pipe resistance should be computed in the process of this call
        self.assertAlmostEqual(bh.calc_bh_resist(flow_rate=0.5, temperature=20), 0.23245, delta=0.0001)
