from math import log as ln
from math import pi

from bhr.enums import DoubleUPipeInletArrangement
from bhr.u_tube import UTube
from bhr.utilities import coth


class DoubleUTube(UTube):

    def __init__(self,
                 borehole_diameter: float,  # m
                 pipe_outer_diameter: float,  # m
                 pipe_dimension_ratio: float,  # unitless, ratio of pipe outer diameter / thickness
                 length: float,  # m
                 shank_space: float,  # distance between adjacent pipe centers, assumes symmetrical placement
                 pipe_conductivity: float,  # W/(m-K)
                 pipe_inlet_arrangement: str,
                 grout_conductivity: float,  # W/(m-K)
                 soil_conductivity: float,  # W/(m-K)
                 fluid_type: str,
                 fluid_concentration: float = 0):

        super().__init__(pipe_outer_diameter, pipe_dimension_ratio, length, shank_space, pipe_conductivity,
                         fluid_type, fluid_concentration)

        # static parameters
        self.grout_conductivity = grout_conductivity
        self.borehole_radius = borehole_diameter / 2  # radius of borehole (m)
        self.pipe_radius = pipe_outer_diameter / 2  # pipe outer radius (m)

        if pipe_inlet_arrangement == DoubleUPipeInletArrangement.ADJACENT.name:
            self.pipe_inlet_arrangement = DoubleUPipeInletArrangement.ADJACENT
        elif pipe_inlet_arrangement == DoubleUPipeInletArrangement.DIAGONAL.name:
            self.pipe_inlet_arrangement = DoubleUPipeInletArrangement.DIAGONAL
        else:
            msg = "Invalid pipe_inlet_arrangement. Use one of the allowed values: " + \
                  f"[{'\", \"'.join(map(str, [DoubleUPipeInletArrangement.ADJACENT.name,
                                              DoubleUPipeInletArrangement.DIAGONAL.name]))}]"
            raise AssertionError(msg)

        self.bh_length = length  # length of borehole (m)
        self.grout_conductivity = grout_conductivity  # W/(m-K)
        self.soil_conductivity = soil_conductivity  # W/(m-K)

        # (m) radial distance between centers of symmetrically placed pipes and borehole center (rc)
        self.pipe_centers_radius = shank_space * (2 ** 0.5 / 2)

        self.sigma = (self.grout_conductivity - self.soil_conductivity) / (
                self.grout_conductivity + self.soil_conductivity)  # thermal conductivity ratio, dimensionless

        # Check if shank spacing realistic
        assert self.pipe_radius * 2 <= shank_space <= 2 / 2 ** 0.5 * (self.borehole_radius - self.pipe_radius), (
                'Shank space is not within bounds. ' +
                'MAX is 2 / sqrt(2) * (borehole_radius - pipe_outer_radius). MIN is pipe_outer_radius * 2')

        # non-static parameters
        self.pipe_resist = None

    def update_b1(self, flow_rate: float, temperature: float) -> float:
        """
        Updates b1 coefficient.

        Javed, S. & Spitler, J.D. Calculation of Borehole Thermal Resistance. In 'Advances in
        Ground-Source Heat Pump Systems,' pp. 84. Rees, S.J. ed. Cambridge, MA. Elsevier Ltd. 2016.

        Eq: 3-47

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Eq: 14

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius
        :return b1: a ratio of (1-beta)/(1+beta), dependent on pipe resistance
                    & grout conductivity, dimensionless
        """

        self.pipe_resist = self.calc_pipe_resist(flow_rate, temperature)
        beta = 2 * pi * self.grout_conductivity * self.pipe_resist
        b1 = (1 - beta) / (1 + beta)  # dimensionless parameter

        return b1

    def calc_bh_resist_local(self, flow_rate: float, temperature: float) -> float:
        """
        Calculates tube-to-borehole resistance (aka local borehole resistance).

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 13 & 14

        :return borehole_resist_local: local borehole resistance K/(W/m)

        """

        b1 = self.update_b1(flow_rate, temperature)

        # static parameters
        p_pc = self.pipe_radius ** 2 / (4 * self.pipe_centers_radius ** 2)  # dimensionless parameter
        p_c = self.pipe_centers_radius ** 2 / (self.borehole_radius ** 8 - self.pipe_centers_radius ** 8) ** (
                1 / 4)  # dimensionless parameter
        p_b = self.borehole_radius ** 2 / (self.borehole_radius ** 8 - self.pipe_centers_radius ** 8) ** (
                1 / 4)  # dimensionless parameter

        # --Borehole resistance, 0th order [K/(W/m)]--
        Rb0 = self.pipe_resist / 4 + 1 / (8 * pi * self.grout_conductivity) * ((ln(self.borehole_radius ** 4 / (
                4 * self.pipe_radius * self.pipe_centers_radius ** 3))) + self.sigma * ln(
            self.borehole_radius ** 8 / (self.borehole_radius ** 8 - self.pipe_centers_radius ** 8)))

        # --Borehole resistance, 1st order [K/(W/m)]--
        borehole_resist_local = Rb0 - 1 / (8 * pi * self.grout_conductivity) * (
                b1 * p_pc * (3 - 8 * self.sigma * p_c ** 4) ** 2
        ) / (1 + b1 * p_pc * (5 + 64 * self.sigma * p_c ** 4 * p_b ** 4))

        return borehole_resist_local

    def calc_internal_resist(self, flow_rate: float, temperature: float) -> float:
        """
        Calculates tube-to-tube resistance (aka internal resistance).

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 18, 19, 22, 23

        :return internal_resist: local internal resistance K/(W/m)
        """

        b1 = self.update_b1(flow_rate, temperature)

        # static parameters
        p_pc = self.pipe_radius ** 2 / (4 * self.pipe_centers_radius ** 2)  # dimensionless parameter
        p_c = self.pipe_centers_radius ** 2 / (self.borehole_radius ** 8 - self.pipe_centers_radius ** 8) ** (
                1 / 4)  # dimensionless parameter
        p_b = self.borehole_radius ** 2 / (self.borehole_radius ** 8 - self.pipe_centers_radius ** 8) ** (
                1 / 4)  # dimensionless parameter

        two_pi_kg = 2 * pi * self.grout_conductivity
        c_1 = self.pipe_centers_radius / self.pipe_radius

        if self.pipe_inlet_arrangement == DoubleUPipeInletArrangement.DIAGONAL:
            # 0th order
            c_2 = self.borehole_radius ** 4 + self.pipe_centers_radius ** 4
            c_3 = self.borehole_radius ** 4 - self.pipe_centers_radius ** 4
            Ra0 = 2 * self.pipe_resist + 2 / two_pi_kg * (ln(c_1) + self.sigma * ln(c_2 / c_3))

            # 1st order

            internal_resist = Ra0 - 2 / two_pi_kg * (
                    b1 * p_pc * (1 + 8 * self.sigma * p_c ** 2 * p_b ** 2) ** 2
            ) / (1 - b1 * p_pc * (
                    3 - 32 * self.sigma * (p_c ** 2 * p_b ** 6 + p_c ** 6 * p_b ** 2)))

            return internal_resist

        elif self.pipe_inlet_arrangement == DoubleUPipeInletArrangement.ADJACENT:
            # 0th order
            d_2 = self.borehole_radius ** 2 + self.pipe_centers_radius ** 2
            d_3 = self.borehole_radius ** 2 - self.pipe_centers_radius ** 2
            Ra0 = 2 * self.pipe_resist + 2 / two_pi_kg * (ln(2 * c_1) + self.sigma * ln(d_2 / d_3))

            # 1st order
            matrix_element_11 = 1 + 16 * b1 * self.sigma * p_pc * (
                    3 * p_c ** 3 * p_b ** 5 + p_c ** 7 * p_b)  # matrix variable
            matrix_element_22 = -1 - 16 * b1 * self.sigma * p_pc * (
                    p_c * p_b ** 7 + 3 * p_c ** 5 * p_b ** 3)  # matrix variable
            matrix_element_21 = b1 * p_pc  # matrix variable
            vector_1 = 1 - 8 * self.sigma * p_c ** 3 * p_b  # vector variable
            vector_2 = 3 + 8 * self.sigma * p_c * p_b ** 3  # vector variable

            internal_resist = Ra0 + 2 / two_pi_kg * b1 * p_pc / 2 * (
                    vector_2 ** 2 * matrix_element_11 - 2 * vector_1 * vector_2 * matrix_element_21 -
                    vector_1 ** 2 * matrix_element_22) / (
                                      matrix_element_11 * matrix_element_22 + matrix_element_21 ** 2)

            return internal_resist

        raise AssertionError("Developer error. Invalid pipe inlet arrangement.")

    def calc_effective_bh_resistance_uhf(self, flow_rate: float, temperature: float) -> float:
        """
        Calculates effective borehole resistance for uniform heat flux along the borehole.

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 44

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius
        :return effective_bhr_uhf: effective borehole resistance under uniform heat flux boundary conditions
                                   [K/(W/m)]
        """

        internal_resist = self.calc_internal_resist(flow_rate, temperature)
        borehole_resist_local = self.calc_bh_resist_local(flow_rate, temperature)

        Rv = self.bh_length / (self.fluid.cp(temperature) * flow_rate)  # (K/(w/m)) thermal resistance factor

        effective_bhr_uhf = borehole_resist_local + Rv ** 2 / (6 * internal_resist)

        return effective_bhr_uhf

    def calc_effective_bh_resistance_ubwt(self, flow_rate: float, temperature: float) -> float:
        """
        Calculates effective borehole resistance for uniform borehole wall temperature.

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 46

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius
        :return effective_bhr_ubwt: effective borehole resistance for uniform borehole wall temperature
                                    boundary condition [K/(W/m)]
        """

        internal_resist = self.calc_internal_resist(flow_rate, temperature)
        borehole_resist_local = self.calc_bh_resist_local(flow_rate, temperature)

        Rv = self.bh_length / (flow_rate * self.fluid.cp(temperature))  # (K/(w/m)) thermal resistance factor
        n = Rv / (2 * borehole_resist_local * internal_resist) ** (1 / 2)
        effective_bhr_ubwt = borehole_resist_local * n * coth(n)

        return effective_bhr_ubwt
