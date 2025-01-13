from unittest import TestCase

from bhr.double_u_borehole import DoubleUTube


class TestDoubleUBorehole(TestCase):
    def setUp(self):
        self.inputs = {
            "borehole_diameter": 0.115,
            "pipe_outer_diameter": 0.032,
            "pipe_dimension_ratio": 11,
            "length": 200,
            "shank_space": 0.032,
            "pipe_conductivity": 0.389,
            "pipe_inlet_arrangement": "DIAGONAL",
            "grout_conductivity": 1.5,
            "soil_conductivity": 3,
            "fluid_type": "WATER"
        }

    def test_init(self):
        bh = DoubleUTube(**self.inputs)
        self.assertEqual(bh.bh_length, 200)

    def test_shank_space_assert(self):
        #checks to make sure code is raising error when shank space is too large
        self.inputs.update({"shank_space": 0.115})
        self.assertRaises(AssertionError) #passes when AssertionError is raised

        #checks to make sure code is raising error when shank space is too small
        self.inputs.update({"shank_space": 0.031})
        self.assertRaises(AssertionError) #passes when AssertionError is raised


    def test_calc_resistances_Diagonal(self):
        bh = DoubleUTube(**self.inputs)

        tolerance = 1e-3
        flow_rate = 0.2077
        temperature = 20
        pipe_resist = 0.05
        # the following tests rely on each previous test
        self.assertAlmostEqual(bh.update_beta_b1(flow_rate, temperature, pipe_resist), 0.359, delta=tolerance)
        self.assertAlmostEqual(bh.calc_bh_resist_local(flow_rate, temperature, pipe_resist), 7.509E-02, delta=tolerance)
        self.assertAlmostEqual(bh.calc_internal_resist_pipe(flow_rate, temperature, pipe_resist), 0.1604,
                               delta=tolerance)
        self.assertAlmostEqual(bh.calc_effective_bh_resistance_uhf(flow_rate, temperature), 0.1302, delta=tolerance)
        self.assertAlmostEqual(bh.calc_effective_bh_resistance_ubwt(flow_rate, temperature), 0.1235, delta=tolerance)


    def test_calc_resistances_Adjacent(self):
        self.inputs.update({"pipe_inlet_arrangement": "ADJACENT"})
        bh = DoubleUTube(**self.inputs)

        tolerance = 1e-3
        flow_rate = 0.2077
        temperature = 20
        pipe_resist = 0.05
        # the following tests rely on each previous test
        self.assertAlmostEqual(bh.update_beta_b1(flow_rate, temperature, pipe_resist), 0.359, delta=tolerance)
        self.assertAlmostEqual(bh.calc_bh_resist_local(flow_rate, temperature, pipe_resist), 7.509E-02, delta=tolerance)
        self.assertAlmostEqual(bh.calc_internal_resist_pipe(flow_rate, temperature, pipe_resist), 0.2617,
                               delta=tolerance)
        self.assertAlmostEqual(bh.calc_effective_bh_resistance_uhf(flow_rate, temperature), 0.1089, delta=tolerance)
        self.assertAlmostEqual(bh.calc_effective_bh_resistance_ubwt(flow_rate, temperature), 0.1062, delta=tolerance)
