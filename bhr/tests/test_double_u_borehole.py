from unittest import TestCase

from bhr.double_u_borehole import DoubleUTube


class TestSingleUBorehole(TestCase):

    def setUp(self):
        self.inputs = {
            "borehole_diameter": 0.096,
            "pipe_outer_diameter": 0.032,
            "pipe_dimension_ratio": 11,
            "pipe_config": "ADJACENT",
            "length": 100,
            "shank_space": 0.032,
            "pipe_conductivity": 0.389,
            "grout_conductivity": 0.6,
            "soil_conductivity": 4.0,
            "fluid_type": "WATER"
        }

    def test_init(self):
        bh = DoubleUTube(**self.inputs)
        self.assertEqual(bh.length, 100)

    # def test_calc_internal_and_grout_resistance(self):
    #     bh = SingleUBorehole(**self.inputs)

    #     tolerance = 1e-3
    #     self.assertAlmostEqual(bh.theta_1, 0.33333, delta=tolerance)
    #     self.assertAlmostEqual(bh.theta_2, 3.0, delta=tolerance)
    #     self.assertAlmostEqual(bh.calc_bh_total_internal_resistance(pipe_resist=0.05), 0.32365, delta=tolerance)
    #     self.assertAlmostEqual(bh.calc_bh_grout_resistance(pipe_resist=0.05), 0.17701, delta=tolerance)

    #     self.inputs.update({"soil_conductivity": 1.0, "grout_conductivity": 3.6})
    #     bh = SingleUBorehole(**self.inputs)
    #     self.assertAlmostEqual(bh.calc_bh_total_internal_resistance(pipe_resist=0.05), 0.17456,
    #                            delta=tolerance)
    #     self.assertAlmostEqual(bh.calc_bh_grout_resistance(pipe_resist=0.05), 0.03373, delta=tolerance)

    def test_calc_bh_resit(self):
            bh = DoubleUTube(**self.inputs)
            tolerance = 1e-3
            self.assertAlmostEqual(bh.calc_bh_resist(0.5, 20), 0.22615, delta=tolerance)
