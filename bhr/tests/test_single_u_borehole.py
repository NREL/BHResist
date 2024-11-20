from unittest import TestCase

from bhr.single_u_borehole import SingleUBorehole


class TestSingleUBorehole(TestCase):

    def test_init(self):
        inputs = {
            "borehole_diameter": 0.096,
            "pipe_outer_diameter": 0.032,
            "pipe_dimension_ratio": 11,
            "length": 100,
            "shank_spacing": 0.032,
            "pipe_conductivity": 0.389,
            "grout_conductivity": 0.6,
            "soil_conductivity": 4.0,
            "fluid_type": "WATER"
        }

        bh = SingleUBorehole(**inputs)

        self.assertEqual(bh.length, 100)

        tolerance = 1e-3
        self.assertAlmostEqual(bh.theta_1, 0.33333, delta=tolerance)
        self.assertAlmostEqual(bh.theta_2, 3.0, delta=tolerance)
        self.assertAlmostEqual(bh.calc_bh_total_internal_resistance(pipe_resist=0.05), 0.32365, delta=tolerance)
