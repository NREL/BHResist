from bhr.pipe import Pipe
from math import pi, log
from bhr.fluid import get_fluid
from bhr.utilities import smoothing_function

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

        self.outer_pipe = Pipe(outer_pipe_outer_diameter, outer_pipe_dimension_ratio, length, outer_pipe_conductivity,
                               fluid_type, fluid_concentration)
        self.inner_pipe = Pipe(inner_pipe_outer_diameter, inner_pipe_dimension_ratio, length, inner_pipe_conductivity,
                               fluid_type, fluid_concentration)

        self.annular_hydraulic_diameter = self.outer_pipe.pipe_inner_diameter - self.inner_pipe.pipe_outer_diameter

    def re_annulus(self, flow_rate, temp):
        """
        Reynolds number for annulus flow

        :param flow_rate: mass flow rate, kg/s
        :param temp: temperature, C
        :return: Reynolds number
        """
        re = 4 * flow_rate / (self.fluid.mu(temp) * pi * self.annular_hydraulic_diameter)
        # print("mu = ", self.fluid.mu(temp))
        # print("annular_hydraulic_diameter = ", self.annular_hydraulic_diameter)
        return re

    def laminar_nusselt_annulus(self):
        """
        Laminar Nusselt numbers for annulus flow

        Hellström, G. 1991. Ground Heat Storage: Thermal Analyses of Duct Storage Systems.
        Department of Mathematical Physics, University of Lund, Sweden., p67-71

        :return: Nusselt number for inner surface of annulus pipe, Nusselt number for outer annulus pipe surface
        """
        Nu_ii = 3.66 + 1.2 * (self.inner_pipe.pipe_outer_diameter/self.outer_pipe.pipe_inner_diameter) ** -0.8
        Nu_oo = 3.66 + 1.2 * (self.inner_pipe.pipe_outer_diameter/self.outer_pipe.pipe_inner_diameter) ** 0.5

        return Nu_ii, Nu_oo

    def turbulent_nusselt_annulus(self, re, temp):
        """
        Turbulent Nusselt numbers for annulus flow

        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Eqns 4.10 and 4.11 based on the Dittus-Boelter equation

        :param flow_rate: mass flow rate, kg/s
        :param temp: temperature, C
        :return: Nusselt number for inner surface of annulus pipe, Nusselt number for outer annulus pipe surface
        """

        pr = self.fluid.prandtl(temp)

        Nu_ii = 0.023 * re ** 0.8 * pr ** 0.35
        Nu_oo = Nu_ii

        return Nu_ii, Nu_oo



    def convective_heat_transfer_coefficients_annulus(self, flow_rate, temp):
        """
        Convective heat transfer coefficients for annulus flow

        :param flow_rate: mass flow rate, kg/s
        :param temp: temperature, C
        :return: annulus pipe convective heat transfer coeffiecients for inner and outer surfaces, W/(m^2K)
        """
        low_reynolds = 2300 #limit determined from Hellstrom, G. 1991. Ground Heat Storage: Thermal Analyses of Duct Storage Systems. Department of Mathmatical Physics, University of Lund, Sweden.
        high_reynolds = 10000 #limit based on dittus-boelter equation

        re = self.re_annulus(flow_rate, temp)

        if re < low_reynolds:
           # use this nusslet number when the flow is laminar
            nu_ii = self.laminar_nusselt_annulus()[0]
            nu_oo = self.laminar_nusselt_annulus()[1]
            print("flow is laminar, Re = ", re)
            # print( "Nu_ii = %f, Nu_oo = %f" % (nu_ii, nu_oo))

        elif low_reynolds <= re < high_reynolds:

            #in between
            nu_ii_low = self.laminar_nusselt_annulus()[0]
            nu_ii_high = self.turbulent_nusselt_annulus(10000, temp)[0]
            sigma = smoothing_function(re, a=6150, b=600)
            nu_ii = (1 - sigma) * nu_ii_low + sigma * nu_ii_high

            nu_oo_low = self.laminar_nusselt_annulus()[1]
            nu_oo_high = self.turbulent_nusselt_annulus(10000, temp)[1]
            sigma = smoothing_function(re, a=6150, b=600)
            nu_oo = (1 - sigma) * nu_oo_low + sigma * nu_oo_high

            print("flow is transitional, Re = ", re)
            # print("Nu_ii_low = ", nu_ii_low)
            # print("Nu_ii_high = ", nu_ii_high)
            # print("sigma = ", sigma)
            # print( "Nu_ii = %f, Nu_oo = %f" % (nu_ii, nu_oo))

        else:
            #use this nusslet number when the flow is fully turbulent
            nu_ii = self.turbulent_nusselt_annulus(re, temp)[0]
            nu_oo = self.turbulent_nusselt_annulus(re, temp)[1]
            print("flow is turbulant, Re = ", re)
            # print( "Nu_ii = %f, Nu_oo = %f" % (nu_ii, nu_oo))

        h_outside_inner_pipe = nu_ii * self.fluid.k(temp) / (self.annular_hydraulic_diameter)
        h_inside_outer_pipe = nu_oo * self.fluid.k(temp) / (self.annular_hydraulic_diameter)
        # print( "h_outside_inner_pipe = ", h_outside_inner_pipe)
        # print( "h_inside_outer_pipe = ", h_inside_outer_pipe)
        return h_outside_inner_pipe, h_inside_outer_pipe

    def calc_bh_resist(self, flow_rate, temp):

        #resistances progressing from inside to outside
        r_conv_inner_pipe = self.inner_pipe.calc_pipe_internal_conv_resist(flow_rate, temp)
        r_cond_inner_pipe = self.inner_pipe.calc_pipe_cond_resist()
        r_conv_outer_pipe_inner_wall = 1/(self.convective_heat_transfer_coefficients_annulus(flow_rate, temp)[0] * self.inner_pipe.pipe_outer_diameter *pi)
        r_conv_outer_pipe_outer_wall = 1/(self.convective_heat_transfer_coefficients_annulus(flow_rate, temp)[1] * self.outer_pipe.pipe_inner_diameter *pi)
        r_cond_outer_pipe = self.outer_pipe.calc_pipe_cond_resist()
        r_cond_grout = log( self.borehole_diameter/self.outer_pipe.pipe_outer_diameter) / (2 * pi * self.grout_conductivity)

        print( "r_conv_inner_pipe = ", r_conv_inner_pipe)
        print( "r_cond_inner_pipe = ", r_cond_inner_pipe)
        print( "r_conv_outer_pipe_inner_wall = ", r_conv_outer_pipe_inner_wall)
        print( "r_conv_outer_pipe_outer_wall = ", r_conv_outer_pipe_outer_wall)
        print( "r_cond_outer_pipe = ", r_cond_outer_pipe)
        print( "r_cond_grout = ", r_cond_grout)

        bh_resist = sum([r_conv_inner_pipe, r_cond_inner_pipe, r_conv_outer_pipe_inner_wall, r_conv_outer_pipe_outer_wall,
                         r_cond_outer_pipe, r_cond_grout])

        return bh_resist