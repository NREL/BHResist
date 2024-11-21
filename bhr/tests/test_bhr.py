import unittest

from bhr.borehole import Borehole


class TestBHR(unittest.TestCase):

    def test_init_single_u_tube_bh(self):
        inputs = {
            "fluid_type": "PROPYLENEGLYCOL",
            "fluid_concentration": 0.2,
            "boundary_condition": "UHF",
            "pipes": {
                "pipe-type": "single-u-tube",
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
                    "shank_space": 0.01,
                    "hydraulic_configuration": "PARALLEL"  # or SERIES
                },
                "coaxial": {
                    "inner_pipe": {
                        "outer_diameter": 0.042,
                        "dimension_ratio": 11,
                        "conductivity": 0.4
                    },
                    "outer-pipe": {
                        "outer_diameter": 0.096,
                        "dimension_ratio": 11,
                        "conductivity": 0.4
                    }
                }
            },
            "grout_conductivity": 1.2,
            "soil_conductivity": 2.5,
            "length": 100,
            "borehole_diameter": 0.14
        }

        bh = Borehole()
        bh.init_from_dict(inputs)

        # only pass flow rate, so pipe resistance should be computed in the process of this call

        self.assertEqual(bh.calc_bh_resist(temperature=20, flow_rate=0.5), 0.5)
