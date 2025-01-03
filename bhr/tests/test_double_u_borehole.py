from unittest import TestCase

from bhr.double_u_borehole import DoubleUTube


class TestDoubleUBorehole(TestCase):

    def setUp(self):
        self.inputs = {
            "borehole_diameter": 0.0575,#meters
            "pipe_outer_diameter": 0.016,#m
            "pipe_dimension_ratio": 11,# ratio of pipe outer diameter/thickness, unitless
            "pipe_config": "DIAGONAL",
            "length": 400,#m length of 1 leg of pipe
            "shank_space": 0.02263, #m
            "pipe_conductivity": 0.389,
            "grout_conductivity": 1.5, #W/(m-K)
            "soil_conductivity": 3.0, #W/(m-K)
            "fluid_type": "WATER"
        }

    def test_init(self):
        bh = DoubleUTube(**self.inputs)
        self.assertEqual(bh.length, 200)

    def test_calc_internal_and_grout_resistance(self):
         bh = DoubleUTube(**self.inputs)

         tolerance = 1e-3
         self.assertAlmostEqual(bh.calc_bh_resist, 0.33333, delta=tolerance)
         self.assertAlmostEqual(bh.calc_internal_resist_pipe, 3.0, delta=tolerance)
         self.assertAlmostEqual(bh.calc_effective_bh_resist_uhf, 0.32365, delta=tolerance)
         self.assertAlmostEqual(bh.calc_effective_bh_resist_ubwt, 0.17701, delta=tolerance)
         self.assertAlmostEqual(bh.calc_effective_bh_resist_ave, 0.17701, delta=tolerance)


    #     self.inputs.update({"soil_conductivity": 1.0, "grout_conductivity": 3.6})
    #     bh = SingleUBorehole(**self.inputs)
    #     self.assertAlmostEqual(bh.calc_bh_total_internal_resistance(pipe_resist=0.05), 0.17456,
    #                            delta=tolerance)
    #     self.assertAlmostEqual(bh.calc_bh_grout_resistance(pipe_resist=0.05), 0.03373, delta=tolerance)

    def test_calc_bh_resit(self):
        bh = DoubleUTube(**self.inputs)
        tolerance = 1e-3
        self.assertAlmostEqual(bh.calc_bh_resist(0.2077, 20), 0.03457, delta=tolerance)
