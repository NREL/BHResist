from math import log, pi

from bhr.fluid import get_fluid
from bhr.utilities import inch_to_m, smoothing_function


class Pipe:

    def __init__(self, inputs: dict):

        self.fluid = get_fluid(inputs.get("fluid-type", "WATER"),
                               inputs.get("fluid-concentration", 0))

        # ratio of outer diameter to wall thickness
        self.dimension_ratio = inputs.get("dimension-ratio", 11)

        # set diameters
        if "nominal-pipe-diameter-inch" in inputs:
            self.inner_diameter, self.outer_diameter = \
                self.get_pipe_diameters_imperial(inputs["nominal-pipe-diameter-inch"],
                                                 self.dimension_ratio)
        elif "actual-pipe-outer-diameter-meter" in inputs:
            self.outer_diameter = inputs["actual-pipe-outer-diameter-meter"]
            self.inner_diameter = self.get_inner_dia(self.outer_diameter, self.dimension_ratio)

        self.length = inputs['length']

        self.conductivity = inputs.get("conductivity", 0.4)

        # compute other dimensions
        self.inner_radius = self.inner_diameter / 2
        self.outer_radius = self.outer_diameter / 2
        self.wall_thickness = self.outer_radius - self.inner_radius

        # compute cross-sectional areas
        self.area_cr_inner = pi / 4 * self.inner_diameter ** 2
        self.area_cr_outer = pi / 4 * self.outer_diameter ** 2
        self.area_cr_pipe = self.area_cr_outer - self.area_cr_inner

        # compute surface areas
        self.area_s_inner = pi * self.inner_diameter * self.length
        self.area_s_outer = pi * self.outer_diameter * self.length

        # compute volumes
        self.total_vol = self.area_cr_outer * self.length
        self.fluid_vol = self.area_cr_inner * self.length
        self.pipe_wall_vol = self.area_cr_pipe * self.length

    @staticmethod
    def get_inner_dia(outer_dia: float, dimension_ratio: float) -> float:
        return outer_dia * (1 - 2 / dimension_ratio)

    def get_pipe_diameters_imperial(self, nominimal_pipe_size_inches: float, dimension_ratio: float):
        if nominimal_pipe_size_inches == 0.75:
            outer_dia = 1.05
        elif nominimal_pipe_size_inches == 1.0:
            outer_dia = 1.315
        elif nominimal_pipe_size_inches == 1.25:
            outer_dia = 1.66
        elif nominimal_pipe_size_inches == 1.5:
            outer_dia = 1.9
        elif nominimal_pipe_size_inches == 2.0:
            outer_dia = 2.375
        elif nominimal_pipe_size_inches == 3.0:
            outer_dia = 3.5
        elif nominimal_pipe_size_inches == 4.0:
            outer_dia = 4.5
        elif nominimal_pipe_size_inches == 6.0:
            outer_dia = 6.625
        elif nominimal_pipe_size_inches == 8.0:
            outer_dia = 8.625
        else:
            raise ValueError("Unsupported pipe size")

        return inch_to_m(self.get_inner_dia(outer_dia, dimension_ratio)), inch_to_m(outer_dia)

    def mdot_to_vdot(self, m_dot: float, temp: float) -> float:
        """
        Computes volumetric flow rate based on mass flow rate.

        :param m_dot: mass flow rate, in kg/s
        :param temp: temperature, in C
        :return: volumetric flow rate, in m3/s
        """

        return m_dot / self.fluid.density(temp)

    def mdot_to_re(self, m_dot: float, temp: float) -> float:
        """
        Computes Reynolds number based on mass flow rate.

        :param m_dot: mass flow rate, in kg/s
        :param temp: temperature, in C
        :return: Reynolds number, dimensionless
        """

        return 4 * m_dot / (self.fluid.mu(temp) * pi * self.inner_diameter)

    def mdot_to_velocity(self, m_dot: float, temp: float) -> float:
        """
        Computes velocity based on mass flow rate.

        :param m_dot: mass flow rate in, kg/s
        :param temp: temperature, in C
        :return: velocity, in m/s
        """

        return m_dot / (self.area_cr_inner * self.fluid.density(temp))

    def friction_factor(self, m_dot: float, temp: float) -> float:
        """
        Calculates the friction factor in smooth tubes

        Petukhov, B.S. 1970. 'Heat transfer and friction in turbulent pipe flow with variable physical properties.'
        In Advances in Heat Transfer, ed. T.F. Irvine and J.P. Hartnett, Vol. 6. New York Academic Press.
        """

        re = self.mdot_to_re(m_dot, temp)

        # limits picked be within about 1% of actual values
        low_reynolds = 1500
        high_reynolds = 5000

        if re < low_reynolds:
            return self.laminar_friction_factor(re)
        elif low_reynolds <= re < high_reynolds:
            # pure laminar flow
            f_low = self.laminar_friction_factor(re)

            # pure turbulent flow
            f_high = self.turbulent_friction_factor(re)
            sigma = smoothing_function(re, a=3000, b=450)
            return (1 - sigma) * f_low + sigma * f_high
        else:
            return self.turbulent_friction_factor(re)

    @staticmethod
    def laminar_friction_factor(re: float):
        """
        Laminar friction factor

        :param re: Reynolds number
        :return: friction factor
        """

        return 64.0 / re

    @staticmethod
    def turbulent_friction_factor(re: float):
        """
        Turbulent friction factor

        Petukhov, B. S. (1970). Advances in Heat Transfer, volume 6, chapter Heat transfer and
        friction in turbulent pipe flow with variable physical properties, pages 503-564.
        Academic Press, Inc., New York, NY.

        :param re: Reynolds number
        :return: friction factor
        """

        return (0.79 * log(re) - 1.64) ** (-2.0)

    def pressure_loss(self, mdot: float, temp: float) -> float:

        """
        Pressure loss in straight pipe

        :param m_dot: mass flow rate, kg/s
        :return: pressure loss, Pa
        """

        if mdot <= 0:
            return 0

        term_1 = self.friction_factor(mdot, temp) * self.length / self.inner_diameter
        term_2 = (self.fluid.density(temp) * self.mdot_to_velocity(mdot, temp) ** 2) / 2

        return term_1 * term_2

    @staticmethod
    def laminar_nusselt():
        """
        Laminar Nusselt number for smooth pipes

        mean(4.36, 3.66)
        :return: Nusselt number
        """
        return 4.01

    def turbulent_nusselt(self, re: float, temp: float):
        """
        Turbulent Nusselt number for smooth pipes

        Gnielinski, V. 1976. 'New equations for heat and mass transfer in turbulent pipe and channel flow.'
        International Chemical Engineering 16(1976), pp. 359-368.

        :param re: Reynolds number
        :param temperature: Temperature, C
        :return: Nusselt number
        """

        f = self.friction_factor(re, temp)
        pr = self.fluid.pr(temp)
        return (f / 8) * (re - 1000) * pr / (1 + 12.7 * (f / 8) ** 0.5 * (pr ** (2 / 3) - 1))

    def calc_cond_resist(self):
        """
        Calculates the pipe radial conduction thermal resistance, in [K/(W/m)].

        Javed, S. and Spitler, J.D. 2017. 'Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers.' Applied Energy. 187: 790-806.

        :return: conduction resistance, K/(W/m)
        """

        return log(self.outer_diameter / self.inner_diameter) / (2 * pi * self.conductivity)

    def calc_conv_resist(self, m_dot: float, temp: float):
        """
        Calculates the pipe internal convection thermal resistance, in [k/(W/m)]

        Gnielinski, V. 1976. 'New equations for heat and mass transfer in turbulent pipe and channel flow.'
        International Chemical Engineering 16(1976), pp. 359-368.

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, C
        :return convection resistance, K/(W/m)
        """

        low_reynolds = 2000
        high_reynolds = 4000

        re = self.mdot_to_re(m_dot, temp)

        if re < low_reynolds:
            nu = self.laminar_nusselt()
        elif low_reynolds <= re < high_reynolds:
            nu_low = self.laminar_nusselt()
            nu_high = self.turbulent_nusselt(re, temp)
            sigma = smoothing_function(re, a=3000, b=150)
            nu = (1 - sigma) * nu_low + sigma * nu_high
        else:
            nu = self.turbulent_nusselt(re, temp)

        return 1 / (nu * pi * self.fluid.k(temp))

    def calc_resist(self, mdot: float, temp: float):
        """
        Calculates the combined conduction and convection pipe resistance

        Javed, S. and Spitler, J.D. 2017. 'Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers.' Applied Energy. 187: 790-806.

        Equation 3
        """

        return self.calc_conv_resist(mdot, temp) + self.calc_cond_resist()
