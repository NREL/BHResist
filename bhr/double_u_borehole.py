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
                 shank_space: float,  # radial distance between pipe centers, assumes radially symmetrical placement
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
        self.borehole_radius = borehole_diameter / 2 * 1000  # radius of borehole (mm) rb
        self.pipe_radius = pipe_outer_diameter / 2 * 1000  # pipe outer radius (mm) rp
        self.pipe_inlet_arrangement = pipe_inlet_arrangement
        self.bh_length = length  # length of borehole is half the length of one pipe (m)
        self.grout_conductivity = grout_conductivity
        self.soil_conductivity = soil_conductivity  # W/(m-K)
        self.shank_space = shank_space * 1000  # (mm) radial distance between centers of symmetrically placed pipes and borehole center rc
        self.sigma = (self.grout_conductivity - self.soil_conductivity) / (
                self.grout_conductivity + self.soil_conductivity) # thermal conductivity ratio, dimensionless

        # -derived variables-

        # non-static parameters
        self.beta = None
        self.borhole_resist_local = None
        self.internal_resist = None
        self.effective_bhr_UHF = None
        self.effective_bhr_UWT = None
        self.pipe_resist = None
        self.b1 = None

    def update_beta_b1(self,
                  flow_rate: float = None,
                  temperature: float = None,
                  pipe_resist: float = None) -> float:
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

        if pipe_resist:
            self.pipe_resist = pipe_resist
        else:
            self.pipe_resist = self.calc_pipe_resist(flow_rate, temperature)

        self.beta = 2 * pi * self.grout_conductivity * self.pipe_resist

        self.b1 = (1 - self.beta) / (1 + self.beta)  # dimensionless parameter

        return self.b1

    def calc_bh_resist(self, flow_rate, temperature, pipe_resist):
        """
        Calculates tube-to-borehole resistance (aka local borehole resistance) .

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 13 & 14

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius
        :param pipe_resist: pipe conduction and convection resistance, K/(W/m). only used for testing
        """

        self.update_beta(flow_rate, temperature, pipe_resist)

        #static parameters
        p_pc = self.pipe_radius ** 2 / (4 * self.shank_space ** 2)  # dimensionless parameter
        p_c = self.shank_space ** 2 / (self.borehole_radius ** 8 - self.shank_space ** 8) ** (
                1 / 4)  # dimensionless parameter

        # --Borehole resistance, 0th order [K/(W/m)]--
        Rb0 = self.pipe_resist / 4 + 1 / (8 * pi * self.grout_conductivity) * (
                (ln(self.borehole_radius ** 4 / (4 * self.pipe_radius * self.shank_space ** 3))) +
                self.sigma * ln(self.borehole_radius ** 8 / (self.borehole_radius ** 8 - self.shank_space ** 8)))

        # --Borehole resistance, 1st order [K/(W/m)]--
        self.borhole_resist_local = Rb0 - 1 / (8 * pi * self.grout_conductivity) * (
                self.b1 * p_pc * (3 - 8 * self.sigma * p_c ** 4) ** 2
        ) / (1 + self.b1 * p_pc * (5 + 64 * self.sigma * p_c ** 4 * p_b ** 4))

        # debugging statements
        # print("Lambda_b = %.3E" % self.grout_conductivity)
        # print("Lambda_soil = %.3E" % self.soil_conductivity)
        # print("Length = %.3E" % self.bh_length)
        # print("rb = %.3E" % self.borehole_radius)
        # print("rp = %.3E" % self.pipe_radius)
        # print("rc = %.3E" % self.shank_space)
        # print("Rb0 = %.3E" % Rb0)
        # print("Rb1 = %.3E" % self.Rb1)
        # print("Mass flow rate = %.3E" % flow_rate)

        return self.borhole_resist_local

    def calc_internal_resist_pipe(self, flow_rate, temperature, pipe_resist):
        """
        Calculates tube-to-tube resistance (aka internal resistance).

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 18, 19, 22, 23

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius
        :param pipe_resist: pipe conduction and convection resistance, K/(W/m). only used for testing
        """

        self.update_beta_b1(flow_rate, temperature, pipe_resist)

        #static parameters
        p_pc = self.pipe_radius ** 2 / (4 * self.shank_space ** 2)  # dimensionless parameter
        p_c = self.shank_space ** 2 / (self.borehole_radius ** 8 - self.shank_space ** 8) ** (
                1 / 4)  # dimensionless parameter
        p_b = self.borehole_radius ** 2 / (self.borehole_radius ** 8 - self.shank_space ** 8) ** (
                1 / 4)  # dimensionless parameter

        if self.pipe_inlet_arrangement == "DIAGONAL":
            # 0th order
            Ra0 = 2 * self.pipe_resist + 2 / (2 * pi * self.grout_conductivity) * (
                    ln(self.shank_space / self.pipe_radius) +
                    self.sigma * ln((self.borehole_radius ** 4 + self.shank_space ** 4) / (
                    self.borehole_radius ** 4 - self.shank_space ** 4)))

            # 1st order
            self.internal_resist = Ra0 - 2 / (2 * pi * self.grout_conductivity) * (
                    self.b1 * p_pc * (1 + 8 * self.sigma * p_c ** 2 * p_b ** 2) ** 2
            ) / (1 - self.b1 * p_pc * (
                    3 - 32 * self.sigma * (p_c ** 2 * p_b ** 6 + p_c ** 6 * p_b ** 2)))

            print("Rad0 = %.3E" % Ra0)
            print("Rad1 = %.3E" % self.internal_resist)

            return self.internal_resist

        elif self.pipe_inlet_arrangement == "ADJACENT":
            # 0th order
            Ra0 = 2 * self.pipe_resist + 2 / (2 * pi * self.grout_conductivity) * (
                    ln(2 * self.shank_space / self.pipe_radius) +
                    self.sigma * ln((self.borehole_radius ** 2 + self.shank_space ** 2) / (
                    self.borehole_radius ** 2 - self.shank_space ** 2)))

            # 1st order
            matrix_element_11 = 1 + 16 * self.b1 * self.sigma * p_pc * (
                    3 * p_c ** 3 * p_b ** 5 + p_c ** 7 * p_b)  # matrix variable
            matrix_element_22 = -1 - 16 * self.b1 * self.sigma * p_pc * (
                    p_c * p_b ** 7 + 3 * p_c ** 5 * p_b ** 3)  # matrix variable
            matrix_element_21 = self.b1 * p_pc  # matrix variable
            vector_1 = 1 - 8 * self.sigma * p_c ** 3 * p_b  # vector variable
            vector_2 = 3 + 8 * self.sigma * p_c * p_b ** 3  # vector variable

            self.internal_resist = Ra0 + 2 / (2 * pi * self.grout_conductivity) * self.b1 * p_pc / 2 * (
                    vector_2 ** 2 * matrix_element_11 - 2 * vector_1 * vector_2 * matrix_element_21 -
                    vector_1 ** 2 * matrix_element_22) / (matrix_element_11 * matrix_element_22 + matrix_element_21 ** 2)

            # debugging print statements
            print("Rad0 = %.3E" % Ra0)
            print("Rad1 = %.3E" % self.internal_resist)

            return self.internal_resist
        else:
            assert False

    def calc_effective_bh_resistance_uhf(self, flow_rate, temperature):
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

        Rv = self.bh_length / (self.fluid.cp(temperature) * flow_rate)  # (K/(w/m)) thermal resistance factor

        self.effective_bhr_UHF = self.borhole_resist_local + Rv ** 2 / (6 * self.internal_resist)

        print("Rv = %.3f" % Rv)
        print("Rb_eff_d_UHF = %.3E" % self.effective_bhr_UHF)

        return self.effective_bhr_UHF

    def calc_effective_bh_resistance_ubwt(self, flow_rate, temperature):
        """
        Calculates effective borehole resistance for uniform borehole wall temperature.

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        Eq: 46

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius
        # :param pipe_resist: pipe conduction and convection resistance, K/(W/m). only used for testing
        """

        Rv = self.bh_length / (flow_rate * self.fluid.cp(temperature))  # (K/(w/m)) thermal resistance factor
        n = Rv / (2 * self.borhole_resist_local * self.internal_resist) ** (1 / 2)
        self.effective_bhr_UWT = self.borhole_resist_local * n * coth(n)

        # print("Ra_eff_UWT = %.3E\n" % self.Rb_eff_UWT)

        return self.effective_bhr_UWT

    def calc_effective_bh_resistance_ave(self):
        """
        Averages effective borehole resistance results for the two boundary conditions
        Uniform Heat Flux and Uniform Borehole Wall temperature

        No arguments??
        """

        self.eff_bhr_ave = (self.effective_bhr_UHF + self.effective_bhr_UWT) / 2.0
        return self.eff_bhr_ave
