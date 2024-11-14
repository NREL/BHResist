import unittest

from bhr.borehole import Borehole


class TestBHR(unittest.TestCase):

    def test_init_single_u_tube_bh(self):

        inputs = {
            "borehole-type": "UTUBE", # or "COAXIAL"
            "num-u-tubes": 1,
            "hydraulic-configuration": "PARALLEL", # SERIES
            "shank-spacing"
            "circulating-fluid": {
                "type": "PROPYLENEGLYCOL", # or "WATER" or "ETHYLALCOHOL" or "METHYLEALCOHOL
                "concentration": 0.2,
            },
            "pipe": {
                "nominal-pipe-diameter-inch": 1.0,
                "dimension-ratio": 11,
                "conductivity": 0.4
            },
            "filling": {
                "type": "GROUT",
                "conductivity": 1.2,
                "heat-capacity": 1600000
            },
            "depth": 100,
        }

        bh = Borehole(inputs)
