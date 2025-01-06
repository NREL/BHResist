from unittest import TestCase

from bhr.double_u_borehole import DoubleUTube


class TestDoubleUBorehole(TestCase):

    def test_init(self):
        self.inputs = {
            "borehole_diameter": 0.115,#meters
            "pipe_outer_diameter": 0.032,#m
            "pipe_dimension_ratio": 11,# ratio of pipe outer diameter/thickness, unitless
            "pipe_config": "DIAGONAL",
            "length": 200,#m length of 1 leg of pipe
            "shank_space": 0.02263, #m
            "pipe_conductivity": 0.389,
            "grout_conductivity": 1.5, #W/(m-K)
            "soil_conductivity": 3.0, #W/(m-K)
            "fluid_type": "WATER"
        }
        bh = DoubleUTube(**self.inputs)
        self.assertEqual(bh.bh_length, 200)

    def test_calc_internal_and_grout_resistance_D(self):
        self.inputs = {
            "borehole_diameter": 0.115,#meters
            "pipe_outer_diameter": 0.032,#m
            "pipe_dimension_ratio": 11,# ratio of pipe outer diameter/thickness, unitless
            "pipe_config": "DIAGONAL",
            "length": 200,#m length of 1 leg of pipe
            "shank_space": 0.02263, #m
            "pipe_conductivity": 0.389,
            "grout_conductivity": 1.5, #W/(m-K)
            "soil_conductivity": 3.0, #W/(m-K)
            "fluid_type": "WATER"
        }
        bh = DoubleUTube(**self.inputs)

        tolerance = 1e-3
        flow_rate = 0.2077
        temperature = 20
        pipe_resist = 0.05
        #the following tests rely on each previous test
        self.assertAlmostEqual(bh.update_b1(flow_rate, temperature, pipe_resist), 0.359, delta=tolerance)
        self.assertAlmostEqual(bh.calc_bh_resist(flow_rate, temperature,pipe_resist), 7.509E-02, delta=tolerance)
        self.assertAlmostEqual(bh.calc_internal_resist_pipe(flow_rate, temperature,pipe_resist), 0.1604, delta=tolerance)
        self.assertAlmostEqual(bh.calc_effective_bh_resist_uhf(flow_rate, temperature), 0.1302, delta=tolerance)
        self.assertAlmostEqual(bh.calc_effective_bh_resist_ubwt(flow_rate, temperature), 0.1235, delta=tolerance)
        self.assertAlmostEqual(bh.calc_effective_bh_resist_ave(), 0.1269, delta=tolerance)




    def test_calc_internal_and_grout_resistance_A(self):
        self.inputs = {
            "borehole_diameter": 0.115,#meters
            "pipe_outer_diameter": 0.032,#m
            "pipe_dimension_ratio": 11,# ratio of pipe outer diameter/thickness, unitless
            "pipe_config": "ADJACENT",
            "length": 200,#m length of 1 leg of pipe
            "shank_space": 0.02263, #m
            "pipe_conductivity": 0.389,
            "grout_conductivity": 1.5, #W/(m-K)
            "soil_conductivity": 3.0, #W/(m-K)
            "fluid_type": "WATER"}

        bh = DoubleUTube(**self.inputs)

        tolerance = 1e-3
        flow_rate = 0.2077
        temperature = 20
        pipe_resist = 0.05
        #the following tests rely on each previous test
        self.assertAlmostEqual(bh.update_b1(flow_rate, temperature, pipe_resist), 0.359, delta=tolerance)
        self.assertAlmostEqual(bh.calc_bh_resist(flow_rate, temperature,pipe_resist), 7.509E-02, delta=tolerance)
        self.assertAlmostEqual(bh.calc_internal_resist_pipe(flow_rate, temperature,pipe_resist), 0.2617, delta=tolerance)
        self.assertAlmostEqual(bh.calc_effective_bh_resist_uhf(flow_rate, temperature), 0.1089, delta=tolerance)
        self.assertAlmostEqual(bh.calc_effective_bh_resist_ubwt(flow_rate, temperature), 0.1062, delta=tolerance)
        self.assertAlmostEqual(bh.calc_effective_bh_resist_ave(), 0.1075, delta=tolerance)
    #     self.inputs.update({"soil_conductivity": 1.0, "grout_conductivity": 3.6})
    #     bh = SingleUBorehole(**self.inputs)
    #     self.assertAlmostEqual(bh.calc_bh_total_internal_resistance(pipe_resist=0.05), 0.17456,
    #                            delta=tolerance)
    #     self.assertAlmostEqual(bh.calc_bh_grout_resistance(pipe_resist=0.05), 0.03373, delta=tolerance)


