import unittest

from bhr.borehole import Borehole


class TestBorehole(unittest.TestCase):

    def test_init_single_u_tube_bh(self):
        inputs = {
            "fluid_type": "PROPYLENEGLYCOL",
            "fluid_concentration": 0.2,
            "boundary_condition": "uniform_heat_flux",
            "borehole_type": "single_u_tube",
            "single-u-tube": {
                "pipe_outer_diameter": 0.042,
                "pipe_dimension_ratio": 11,
                "pipe_conductivity": 0.4,
                "shank_space": 0.02,
            },
            "double-u-tube": {
                "pipe_outer_diameter": 0.042,
                "pipe_dimension_ratio": 11,
                "pipe_conductivity": 0.4,
                "shank_space": 0.02,
                "hydraulic_configuration": "PARALLEL"  # or SERIES
            },
            "coaxial": {
                "outer_pipe_outer_diameter": 0.042,
                "outer_pipe_dimension_ratio": 11,
                "outer_pipe_conductivity": 0.4,
                "inner_pipe_outer_diameter": 0.096,
                "inner_pipe_dimension_ratio": 11,
                "inner_pipe_conductivity": 0.4
            },
            "grout_conductivity": 1.2,
            "soil_conductivity": 2.5,
            "length": 100,
            "borehole_diameter": 0.14
        }

        bh = Borehole()
        bh.init_from_dict(inputs)

        # only pass flow rate, so pipe resistance should be computed in the process of this call

        self.assertAlmostEqual(bh.calc_bh_resist(temperature=20, flow_rate=0.5), 0.20304, delta=0.0001)
