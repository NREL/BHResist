import unittest

from bhr.borehole import Borehole


class TestBHR(unittest.TestCase):

    def test_init_bh_init(self):

        inputs = {
            "borehole-type": "SINGLEUTUBE",
            "fluid": {
                "fluid-type": "PROPYLENEGLYCOL",
                "fluid-concentration": 0.2,
            },
            "pipe": {
                "nominal-pipe-diameter-inch": 1.0,
                "dimension-ratio": 11,
                "conductivity": 0.4
            },
            "filling": {
                "grout-type": "GROUT",
                "conductivity": 1.2,
                "heat-capacity": 1600000
            },
            "depth": 100,
        }

        bh = Borehole(inputs)
