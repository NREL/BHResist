from math import pi, log

from bhr.fluid import get_fluid
from bhr.pipe import Pipe
from bhr.utilities import smoothing_function, coth


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
                 fluid_concentration: float):

        self.borehole_diameter = borehole_diameter
        self.grout_conductivity = grout_conductivity
        self.soil_conductivity = soil_conductivity
        self.fluid = get_fluid(fluid_type, fluid_concentration)
        self.length = length

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

        return 4 * flow_rate / (self.fluid.mu(temp) * pi * self.annular_hydraulic_diameter)

    def laminar_nusselt_annulus(self):
        """
        Laminar Nusselt numbers for annulus flow

        Hellström, G. 1991. Ground Heat Storage: Thermal Analyses of Duct Storage Systems.
        Department of Mathematical Physics, University of Lund, Sweden. pp 67-71

        :return nu_ii: Laminar Nusselt number for inner surface of annulus pipe
        :return nu_oo: Laminar Nusselt number for outer annulus pipe surface
        """
        nu_ii = 3.66 + 1.2 * (self.inner_pipe.pipe_outer_diameter / self.outer_pipe.pipe_inner_diameter) ** -0.8
        nu_oo = 3.66 + 1.2 * (self.inner_pipe.pipe_outer_diameter / self.outer_pipe.pipe_inner_diameter) ** 0.5

        return nu_ii, nu_oo

    def turbulent_nusselt_annulus(self, re, temp):
        """
        Turbulent Nusselt numbers for annulus flow

        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Eqns 4.10 and 4.11 based on the Dittus-Boelter equation

        :param re: Reynolds number
        :param temp: temperature, C
        :return nu_ii: Turbulent Nusselt number for inner surface of annulus pipe
        :return nu_oo: Turbulent Nusselt number for outer annulus pipe surface
        """

        pr = self.fluid.prandtl(temp)

        nu_ii = 0.023 * re ** 0.8 * pr ** 0.35
        nu_oo = nu_ii

        return nu_ii, nu_oo

    def convective_resist_annulus(self, flow_rate, temp):
        """
        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Eqns 4.4 - 4.11

        :param flow_rate: mass flow rate, kg/s
        :param temp: temperature, C
        :return r_conv_outside_inner_pipe: convective resistances along the outer wall of the inner pipe, K/(W/m)
        :return r_conv_inside_outer_pipe: convective resistance along the inside wall of the outer pipe K/(W/m)
        """

        # limit determined from Hellström, G. 1991. Ground Heat Storage: Thermal Analyses of
        # Duct Storage Systems. Department of Mathematical Physics, University of Lund, Sweden.

        low_reynolds = 2300

        # limit based on Dittus-Boelter equation
        high_reynolds = 10000

        re = self.re_annulus(flow_rate, temp)

        if re < low_reynolds:
            # use this Nusslet number when the flow is laminar
            nu_ii, nu_oo = self.laminar_nusselt_annulus()

        elif low_reynolds <= re < high_reynolds:

            # in between
            nu_ii_low, nu_oo_low = self.laminar_nusselt_annulus()
            nu_ii_high, nu_oo_high = self.turbulent_nusselt_annulus(high_reynolds, temp)
            nu_ii = smoothing_function(re, low_reynolds, high_reynolds, nu_ii_low, nu_ii_high)
            nu_oo = smoothing_function(re, low_reynolds, high_reynolds, nu_oo_low, nu_oo_high)

        else:
            # use this Nusslet number when the flow is fully turbulent
            nu_ii, nu_oo = self.turbulent_nusselt_annulus(re, temp)

        r_conv_outside_inner_pipe = self.annular_hydraulic_diameter / (
                nu_ii * self.fluid.k(temp) * self.inner_pipe.pipe_outer_diameter * pi)

        r_conv_inside_outer_pipe = self.annular_hydraulic_diameter / (
                nu_oo * self.fluid.k(temp) * self.outer_pipe.pipe_inner_diameter * pi)

        return r_conv_outside_inner_pipe, r_conv_inside_outer_pipe

    def calc_local_bh_resistance(self, flow_rate, temp):
        """
        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Eqns 4.4 and 4.5

        :param flow_rate: mass flow rate, kg/s
        :param temp: temperature, C
        :return local_bh_resist: total local borehole resistance K /(W/m)
        :return r_internal_resist: local internal borehole resistance K /(W/m)
        :return r_borehole_resist: local borehole resistance K /(W/m)

        """
        # resistances progressing from inside to outside
        r_conv_inner_pipe = self.inner_pipe.calc_pipe_internal_conv_resist(flow_rate, temp)
        r_cond_inner_pipe = self.inner_pipe.calc_pipe_cond_resist()
        r_conv_outside_inner_pipe = self.convective_resist_annulus(flow_rate, temp)[0]
        r_conv_inside_outer_pipe = self.convective_resist_annulus(flow_rate, temp)[1]
        r_cond_outer_pipe = self.outer_pipe.calc_pipe_cond_resist()
        r_cond_grout = log(self.borehole_diameter / self.outer_pipe.pipe_outer_diameter) / (
                2 * pi * self.grout_conductivity)

        r_internal_resist = sum([r_conv_inner_pipe, r_cond_inner_pipe, r_conv_outside_inner_pipe])
        r_borehole_resist = sum([r_conv_inside_outer_pipe,
                                 r_cond_outer_pipe, r_cond_grout])
        local_bh_resist = r_internal_resist + r_borehole_resist

        return [local_bh_resist, r_internal_resist, r_borehole_resist]

    def calc_effective_bh_resistance_uhf(self, flow_rate, temp):
        """
        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Eqn 4.33

        :param flow_rate: mass flow rate, kg/s
        :param temp: temperature, C
        :return effective_bhr_uhf: effective borehole resistance for
                 uniform heat flux boundary condition  [K/(W/m)]
        """

        _, r_a, r_b = self.calc_local_bh_resistance(flow_rate, temp)
        Rv = self.length / (flow_rate * self.fluid.cp(temp))  # (K/(w/m)) thermal resistance factor
        effective_bhr_uhf = r_b + 1 / (3 * r_a) * (Rv) ** 2

        return effective_bhr_uhf

    def calc_effective_bh_resistance_ubwt(self, flow_rate, temp):
        """
        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Eqns 4.28 & 4.29

        :param flow_rate: mass flow rate, kg/s
        :param temp: temperature, C
        :return effective_bhr_ubwt: effective borehole resistance for
                 uniform borehole wall temperature boundary condition  [K/(W/m)]
        """

        _, r_a, r_b = self.calc_local_bh_resistance(flow_rate, temp)
        Rv = self.length / (flow_rate * self.fluid.cp(temp))  # (K/(w/m)) thermal resistance factor
        n = Rv / (2 * r_b) * (1 + 4 * r_b / r_a) ** (1 / 2)
        effective_bhr_ubwt = r_b * n * coth(n)

        return effective_bhr_ubwt
