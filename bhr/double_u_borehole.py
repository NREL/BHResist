from math import pi, log

from bhr.u_tube import UTube
from bhr.utilities import coth

# math functions
ln = log


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
        self.borehole_radius = borehole_diameter / 2  # radius of borehole (m) rb
        self.pipe_radius = pipe_outer_diameter / 2  # pipe outer radius (m) rp
        self.pipe_inlet_arrangement = pipe_inlet_arrangement
        self.bh_length = length  # length of borehole is half the length of one pipe (m)
        self.grout_conductivity = grout_conductivity
        self.soil_conductivity = soil_conductivity  # W/(m-K)
        self.pipe_centers_radius = shank_space * (
                    2 ** 0.5 / 2)  # (m) radial distance between centers of symmetrically placed pipes and borehole center (rc)
        self.sigma = (self.grout_conductivity - self.soil_conductivity) / (
                self.grout_conductivity + self.soil_conductivity)  # thermal conductivity ratio, dimensionless

        # Check if shank spacing realistic
        assert self.pipe_radius * 2 <= shank_space <= 2 / 2 ** 0.5 * (self.borehole_radius - self.pipe_radius), (
                    'Shank space is not within bounds. ' +
                    'MAX is 2 / sqrt(2) * (borehole_radius - pipe_outer_radius). MIN is pipe_outer_radius * 2')

        # non-static parameters
        self.beta = None
        self.effective_bhr_UHF = None
        self.effective_bhr_UWT = None
        self.pipe_resist = None
        self.b1 = None

    def update_beta_b1(self,
                       flow_rate: float,
                       temperature: float) -> float:
        """
        Updates Beta & b1 coefficients.

        Javed, S. & Spitler, J.D. Calculation of Borehole Thermal Resistance. In 'Advances in
        Ground-Source Heat Pump Systems,' pp. 84. Rees, S.J. ed. Cambridge, MA. Elsevier Ltd. 2016.

        Eq: 3-47

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Eq: 14

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius
        :param pipe_resist: pipe conduction and convection resistance, K/(W/m). only used for testing
        """

        self.pipe_resist = self.calc_pipe_resist(flow_rate, temperature)
        self.beta = 2 * pi * self.grout_conductivity * self.pipe_resist
        self.b1 = (1 - self.beta) / (1 + self.beta)  # dimensionless parameter

        return self.b1

    def calc_bh_resist_local(self):
        """
        Calculates tube-to-borehole resistance (aka local borehole resistance) .

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 13 & 14
        """

        # static parameters
        p_pc = self.pipe_radius ** 2 / (4 * self.pipe_centers_radius ** 2)  # dimensionless parameter
        p_c = self.pipe_centers_radius ** 2 / (self.borehole_radius ** 8 - self.pipe_centers_radius ** 8) ** (
                1 / 4)  # dimensionless parameter
        p_b = self.borehole_radius ** 2 / (self.borehole_radius ** 8 - self.pipe_centers_radius ** 8) ** (
                1 / 4)  # dimensionless parameter

        # --Borehole resistance, 0th order [K/(W/m)]--
        Rb0 = self.pipe_resist / 4 + 1 / (8 * pi * self.grout_conductivity) * (
                (ln(self.borehole_radius ** 4 / (4 * self.pipe_radius * self.pipe_centers_radius ** 3))) +
                self.sigma * ln(
            self.borehole_radius ** 8 / (self.borehole_radius ** 8 - self.pipe_centers_radius ** 8)))

        # --Borehole resistance, 1st order [K/(W/m)]--
        borehole_resist_local = Rb0 - 1 / (8 * pi * self.grout_conductivity) * (
                self.b1 * p_pc * (3 - 8 * self.sigma * p_c ** 4) ** 2
        ) / (1 + self.b1 * p_pc * (5 + 64 * self.sigma * p_c ** 4 * p_b ** 4))

        return borehole_resist_local

    def calc_internal_resist(self):
        """
        Calculates tube-to-tube resistance (aka internal resistance).

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 18, 19, 22, 23
        """

        # static parameters
        p_pc = self.pipe_radius ** 2 / (4 * self.pipe_centers_radius ** 2)  # dimensionless parameter
        p_c = self.pipe_centers_radius ** 2 / (self.borehole_radius ** 8 - self.pipe_centers_radius ** 8) ** (
                1 / 4)  # dimensionless parameter
        p_b = self.borehole_radius ** 2 / (self.borehole_radius ** 8 - self.pipe_centers_radius ** 8) ** (
                1 / 4)  # dimensionless parameter

        if self.pipe_inlet_arrangement == "DIAGONAL":
            # 0th order
            Ra0 = 2 * self.pipe_resist + 2 / (2 * pi * self.grout_conductivity) * (
                    ln(self.pipe_centers_radius / self.pipe_radius) +
                    self.sigma * ln((self.borehole_radius ** 4 + self.pipe_centers_radius ** 4) / (
                    self.borehole_radius ** 4 - self.pipe_centers_radius ** 4)))

            # 1st order
            internal_resist = Ra0 - 2 / (2 * pi * self.grout_conductivity) * (
                    self.b1 * p_pc * (1 + 8 * self.sigma * p_c ** 2 * p_b ** 2) ** 2
            ) / (1 - self.b1 * p_pc * (
                    3 - 32 * self.sigma * (p_c ** 2 * p_b ** 6 + p_c ** 6 * p_b ** 2)))

            return internal_resist

        elif self.pipe_inlet_arrangement == "ADJACENT":
            # 0th order
            Ra0 = 2 * self.pipe_resist + 2 / (2 * pi * self.grout_conductivity) * (
                    ln(2 * self.pipe_centers_radius / self.pipe_radius) +
                    self.sigma * ln((self.borehole_radius ** 2 + self.pipe_centers_radius ** 2) / (
                    self.borehole_radius ** 2 - self.pipe_centers_radius ** 2)))

            # 1st order
            matrix_element_11 = 1 + 16 * self.b1 * self.sigma * p_pc * (
                    3 * p_c ** 3 * p_b ** 5 + p_c ** 7 * p_b)  # matrix variable
            matrix_element_22 = -1 - 16 * self.b1 * self.sigma * p_pc * (
                    p_c * p_b ** 7 + 3 * p_c ** 5 * p_b ** 3)  # matrix variable
            matrix_element_21 = self.b1 * p_pc  # matrix variable
            vector_1 = 1 - 8 * self.sigma * p_c ** 3 * p_b  # vector variable
            vector_2 = 3 + 8 * self.sigma * p_c * p_b ** 3  # vector variable

            internal_resist = Ra0 + 2 / (2 * pi * self.grout_conductivity) * self.b1 * p_pc / 2 * (
                    vector_2 ** 2 * matrix_element_11 - 2 * vector_1 * vector_2 * matrix_element_21 -
                    vector_1 ** 2 * matrix_element_22) / (
                                          matrix_element_11 * matrix_element_22 + matrix_element_21 ** 2)

            return internal_resist
        else:
            assert False

    def calc_effective_bh_resistance_uhf(self,
                                         flow_rate,
                                         temperature):
        """
        Calculates effective borehole resistance for uniform heat flux along the borehole.

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 44

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius
        # :param pipe_resist: pipe conduction and convection resistance, K/(W/m). only used for testing
        """

        self.update_beta_b1(flow_rate, temperature)
        internal_resist = self.calc_internal_resist()
        borehole_resist_local = self.calc_bh_resist_local()

        Rv = self.bh_length / (self.fluid.cp(temperature) * flow_rate)  # (K/(w/m)) thermal resistance factor


        self.effective_bhr_UHF = borehole_resist_local + Rv ** 2 / (6 * internal_resist)

        return self.effective_bhr_UHF

    def calc_effective_bh_resistance_ubwt(self,
                                          flow_rate,
                                          temperature):
        """
        Calculates effective borehole resistance for uniform borehole wall temperature.

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 46

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius
        """

        self.update_beta_b1(flow_rate, temperature)

        internal_resist = self.calc_internal_resist()
        borehole_resist_local = self.calc_bh_resist_local()

        Rv = self.bh_length / (flow_rate * self.fluid.cp(temperature))  # (K/(w/m)) thermal resistance factor
        n = Rv / (2 * borehole_resist_local * internal_resist) ** (1 / 2)
        self.effective_bhr_UWT = borehole_resist_local * n * coth(n)

        return self.effective_bhr_UWT
