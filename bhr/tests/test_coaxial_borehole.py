from unittest import TestCase

from bhr.coaxial_borehole import Coaxial


class TestCoaxialBorehole(TestCase):
    def setUp(self):
        self.inputs = {
                 "borehole_diameter": 0.115,
                 "outer_pipe_outer_diameter": 0.064,
                 "outer_pipe_dimension_ratio": 11,
                 "outer_pipe_conductivity": 0.389,
                 "inner_pipe_outer_diameter": 0.032,
                 "inner_pipe_dimension_ratio": 11,
                 "inner_pipe_conductivity": 0.389,
                 "length": 200,
                 "grout_conductivity": 1.5,
                 "soil_conductivity": 3,
                 "fluid_type": "WATER",
                 "fluid_concentration": 0
        }

    def test_init(self):
        coax = Coaxial(**self.inputs)
        self.assertEqual(coax.borehole_diameter, 0.115)

    def test_re_annulus(self):
        coax = Coaxial(**self.inputs)
        self.assertAlmostEqual(coax.re_annulus(flow_rate=0.5,temp=20), 31200.1777, delta=1e-3)

    def test_laminar_nusselt_annulus(self):
        coax = Coaxial(**self.inputs)

        tolerance = 1e-3

        self.assertAlmostEqual(coax.laminar_nusselt_annulus()[0], 5.4394 , delta=tolerance)
        self.assertAlmostEqual(coax.laminar_nusselt_annulus()[1], 4.5980 , delta=tolerance)

    def test_turbulent_nusselt_annulus(self):

         coax = Coaxial(**self.inputs)
         self.assertAlmostEqual(coax.turbulent_nusselt_annulus(flow_rate = 0.5 , temp = 20)[0], 179.02132 , delta=1e-3)
         self.assertAlmostEqual(coax.turbulent_nusselt_annulus(flow_rate = 0.5 , temp = 20)[1], 179.02132 , delta=1e-3)

    def test_convective_heat_transfer_coefficients_annulus(self):
         coax = Coaxial(**self.inputs)

         #tests turbulant flow, convective heat transfer coefficient of outside surface of inner pipe
         self.assertAlmostEqual(coax.convective_heat_transfer_coefficients_annulus(flow_rate = 0.5 , temp = 20)[0], 5260.2484 ,delta=1e-3)
         #tests transitional flow, convective heat transfer coefficient of outside surface of inner pipe
         self.assertAlmostEqual(coax.convective_heat_transfer_coefficients_annulus(flow_rate = 0.1 , temp = 20)[0], 9946168.6213 ,delta=1e-3) #ACK this is a very high value
         #tests laminar flow, convective heat transfer coefficient of outside surface of inner pipe
         self.assertAlmostEqual(coax.convective_heat_transfer_coefficients_annulus(flow_rate = 0.02 , temp = 20)[0], 159.829 ,delta=1e-3)

         #tests turbulant flow, convective heat transfer coefficient of inside of outer pipe
         self.assertAlmostEqual(coax.convective_heat_transfer_coefficients_annulus(flow_rate = 0.5 , temp = 20)[1], 5260.2484 ,delta=1e-3)
         #tests transitional flow, convective heat transfer coefficient of inside of outer pipe
         self.assertAlmostEqual(coax.convective_heat_transfer_coefficients_annulus(flow_rate = 0.1 , temp = 20)[1], 9946168.60286 ,delta=1e-3) #ACK this is a very high value
         #tests laminar flow, convective heat transfer coefficient of inside of outer pipe
         self.assertAlmostEqual(coax.convective_heat_transfer_coefficients_annulus(flow_rate = 0.02 , temp = 20)[1], 135.10714 ,delta=1e-3)

    def test_calc_bh_resist(self):
         coax = Coaxial(**self.inputs)

        #test turbulant result
         self.assertAlmostEqual(coax.calc_bh_resist(flow_rate = 0.5 , temp = 20), 0.23245 , delta=1e-3)
        #test intermediate result
         self.assertAlmostEqual(coax.calc_bh_resist(flow_rate = 0.1 , temp = 20), 0.24006 , delta=1e-3)
         #test laminar result
         self.assertAlmostEqual(coax.calc_bh_resist(flow_rate = 0.02 , temp = 20), 0.4663 , delta=1e-3)