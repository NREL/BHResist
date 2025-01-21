from bhr.pipe import Pipe
from math import pi, log
import bhr.fluid import get_fluid

class Coaxial:

    def __init__(self,
                 borehole_diameter: float,
                 outer_pipe_outer_diameter: float,
                 outer_pipe_dimension_ratio: float,
                 outer_pipe_conductivity: float,
                 inner_pipe_outer_diameter: float,
                 inner_pipe_dimension_ratio: float,
                 inner_pipe_conductivity: float,
                 length: float,
                 grout_conductivity: float,
                 soil_conductivity: float,
                 fluid_type: str,
                 fluid_concentration: float,
                 ):
        self.borehole_diameter = borehole_diameter
        self.grout_conductivity = grout_conductivity
        self.soil_conductivity = soil_conductivity
        self.fluid = get_fluid(fluid_type, fluid_concentration)

        self.inner_pipe = Pipe(outer_pipe_outer_diameter, outer_pipe_dimension_ratio, length, outer_pipe_conductivity,
                               fluid_type, fluid_concentration)
        self.outer_pipe = Pipe(inner_pipe_outer_diameter, inner_pipe_dimension_ratio, length, inner_pipe_conductivity,
                               fluid_type, fluid_concentration)

        self.annular_hydraulic_diameter = self.outer_pipe.pipe_inner_diameter - self.inner_pipe.pipe_outer_diameter

    def convective_heat_transfer_coefficients_concentric_annulus(self, flow_rate, temp):

        low_reynolds = 2000
        high_reynolds = 4000

        re = Pipe.mdot_to_re(flow_rate, temp)
        pr = self.fluid.prandtl(temp)

        if re < low_reynolds:
            Nu_ii = 3.66 + 1.2 (self.inner_pipe.pipe_outer_diameter/self.outer_pipe.pipe_inner_diameter) ** -0.8
            Nu_oo = 3.66 + 1.2 (self.inner_pipe.pipe_inner_diameter/self.outer_pipe.pipe_inner_diameter) ** 0.5
        elif low_reynolds <= re < high_reynolds:
            Nu_ii = 1 #ERROR
            Nu_oo = 1 #ERROR
        else re >= high_reynolds:
            Nu_ii = 0.023 * re ** 0.8 * pr ** 0.35
            Nu_oo = Nu_ii

        h_outside_inner_pipe = Nu_ii * self.fluid.k(temp) / (self.annular_hydraulic_diameter) #hop
        h_inside_annulus = Nu_oo * self.fluid.k(temp) / (self.annular_hydraulic_diameter) #hoa

        return h_outside_inner_pipe, h_inside_annulus

    def calc_bh_resist(self, flow_rate, temp):

        #resistances progressing from inside to outside
        r_conv_inner_pipe = self.inner_pipe.calc_pipe_internal_conv_resist(flow_rate, temp)
        r_cond_inner_pipe = self.inner_pipe.calc_pipe_cond_resist()
        r_conv_outer_pipe_inner_wall = 1/(self.convective_heat_transfer_coefficients_concentric_annulus(flow_rate, temp)[0] * self.inner_pipe.pipe_outer_diameter *pi)
        r_conv_outer_pipe_outer_wall = 1/(self.convective_heat_transfer_coefficients_concentric_annulus(flow_rate, temp)[1] * self.outer_pipe.pipe_inner_diameter *pi)
        r_cond_outer_pipe = self.outer_pipe.calc_pipe_cond_resist()
        r_cond_grout = log(
                self.borehole_diameter/self.outer_pipe.pipe_outer_diameter) / (
                2 * pi * self.grout_conductivity)

        bh_resist = sum([r_conv_inner_pipe, r_cond_inner_pipe, r_conv_outer_pipe_inner_wall, r_conv_outer_pipe_outer_wall,
                         r_cond_outer_pipe, r_cond_grout])

        return bh_resist