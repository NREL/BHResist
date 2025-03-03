from math import log, pi

from bhr.u_tube import UTube
from bhr.utilities import coth


class SingleUBorehole(UTube):

    def __init__(self,
                 borehole_diameter: float,
                 pipe_outer_diameter: float,
                 pipe_dimension_ratio: float,
                 length: float,
                 shank_space: float,
                 pipe_conductivity: float,
                 grout_conductivity: float,
                 soil_conductivity: float,
                 fluid_type: str,
                 fluid_concentration: float = 0):
        super().__init__(pipe_outer_diameter, pipe_dimension_ratio, length, shank_space, pipe_conductivity,
                         fluid_type, fluid_concentration)

        # other static parameters
        self.borehole_diameter = borehole_diameter
        self.grout_conductivity = grout_conductivity
        self.soil_conductivity = soil_conductivity
        self.theta_1 = self.shank_space / self.borehole_diameter
        self.theta_2 = self.borehole_diameter / self.pipe_outer_diameter
        self.theta_3 = 1 / (2 * self.theta_1 * self.theta_2)
        self.sigma = (self.grout_conductivity - self.soil_conductivity) / (
                self.grout_conductivity + self.soil_conductivity)
        self.bh_length = length

        # non-static parameters
        self.pipe_resist = None

    def update_beta(self, flow_rate: float, temperature: float) -> None:
        """
        Updates Beta coefficient.

        Javed, S. & Spitler, J.D. Calculation of Borehole Thermal Resistance. In 'Advances in
        Ground-Source Heat Pump Systems,' pp. 84. Rees, S.J. ed. Cambridge, MA. Elsevier Ltd. 2016.

        Eq: 3-47

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Eq: 14

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius

        :return: none
        """

        self.pipe_resist = self.calc_pipe_resist(flow_rate, temperature)
        beta = 2 * pi * self.grout_conductivity * self.pipe_resist

        return beta

    def calc_direct_coupling_resistance(self, flow_rate, temperature) -> tuple:
        r_a = self.calc_total_internal_bh_resistance(flow_rate, temperature)
        r_b = self.calc_local_bh_resistance(flow_rate, temperature)

        r_12 = (4 * r_a * r_b) / (4 * r_b - r_a)

        # reset if negative
        if r_12 < 0:
            r_12 = 70

        resist_bh_direct_coupling = r_12
        return resist_bh_direct_coupling, r_b

    def calc_local_bh_resistance(self, flow_rate, temperature) -> float:
        """
        Calculates the average thermal resistance of the borehole using the first-order multipole method.

        Resistance between the fluid in the U-tube(s) to the borehole wall (m-K/W)

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Equation 13

        :return: average thermal resistance, K/(W/m)
        """
        beta = self.update_beta(flow_rate, temperature)

        final_term_1 = log(self.theta_2 / (2 * self.theta_1 * (1 - self.theta_1 ** 4) ** self.sigma))

        term_2_num = self.theta_3 ** 2 * (1 - (4 * self.sigma * self.theta_1 ** 4) / (1 - self.theta_1 ** 4)) ** 2
        term_2_den_pt_1 = (1 + beta) / (1 - beta)
        term_2_den_pt_2 = self.theta_3 ** 2 * (1 + (16 * self.sigma * self.theta_1 ** 4) / (1 - self.theta_1 ** 4) ** 2)
        term_2_den = term_2_den_pt_1 + term_2_den_pt_2
        final_term_2 = term_2_num / term_2_den

        resist_bh_ave = (1 / (4 * pi * self.grout_conductivity)) * (beta + final_term_1 - final_term_2)
        return resist_bh_ave

    def calc_total_internal_bh_resistance(self, flow_rate, temperature) -> float:
        """
        Calculates the total internal thermal resistance of the borehole using the first-order multipole method.

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Equation 26

        :return: total internal thermal resistance, K/(W/m)
        """
        beta = self.update_beta(flow_rate, temperature)

        term_1_num = (1 + self.theta_1 ** 2) ** self.sigma
        term_1_den = self.theta_3 * (1 - self.theta_1 ** 2) ** self.sigma
        final_term_1 = log(term_1_num / term_1_den)

        term_2_num = self.theta_3 ** 2 * (1 - self.theta_1 ** 4 + 4 * self.sigma * self.theta_1 ** 2) ** 2
        term_2_den_pt_1 = (1 + beta) / (1 - beta) * (1 - self.theta_1 ** 4) ** 2
        term_2_den_pt_2 = self.theta_3 ** 2 * (1 - self.theta_1 ** 4) ** 2
        term_2_den_pt_3 = 8 * self.sigma * self.theta_1 ** 2 * self.theta_3 ** 2 * (1 + self.theta_1 ** 4)
        term_2_den = term_2_den_pt_1 - term_2_den_pt_2 + term_2_den_pt_3
        final_term_2 = term_2_num / term_2_den

        resist_bh_total_internal = 1 / (pi * self.grout_conductivity) * (beta + final_term_1 - final_term_2)

        return resist_bh_total_internal

    def calc_grout_resistance(self, flow_rate, temperature) -> float:
        """
        Calculates grout resistance. Use for validation.

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Eq: 3

        :return: grout resistance, K/(W-m)
        """

        resist_bh_grout = self.calc_local_bh_resistance(flow_rate,temperature) - self.pipe_resist / 2.0
        return resist_bh_grout

    def calc_effective_bh_resistance_uhf(self, flow_rate: float, temperature: float) -> float:
        """
        Calculates the effective thermal resistance of the borehole assuming a uniform heat flux.

        Javed, S. & Spitler, J.D. Calculation of Borehole Thermal Resistance. In 'Advances in
        Ground-Source Heat Pump Systems,' pp. 84. Rees, S.J. ed. Cambridge, MA. Elsevier Ltd. 2016.

        Eq: 3-67

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius

        :return: effective thermal resistance, K/(W/m)
        """

        r_a = self.calc_total_internal_bh_resistance(flow_rate,temperature)
        r_b = self.calc_local_bh_resistance(flow_rate, temperature)

        pt_1 = 1 / (3 * r_a)
        pt_2 = (self.bh_length / (self.fluid.cp(temperature) * flow_rate)) ** 2
        resist_short_circuiting = pt_1 * pt_2

        resist_bh_effective_uhf = r_b + resist_short_circuiting
        return resist_bh_effective_uhf

    def calc_effective_bh_resistance_ubwt(self, flow_rate: float, temperature: float) -> float:
        """
        Calculates the effective thermal resistance of the borehole assuming a uniform borehole wall temperature.

        Javed, S. & Spitler, J.D. Calculation of Borehole Thermal Resistance. In 'Advances in
        Ground-Source Heat Pump Systems,' pp. 84. Rees, S.J. ed. Cambridge, MA. Elsevier Ltd. 2016.

        :param flow_rate: mass flow rate, kg/s
        :param temperature: temperature, Celsius

        :return: effective thermal resistance, K/(W/m)
        """
        r_a = self.calc_total_internal_bh_resistance(flow_rate,temperature) #R_a
        r_b = self.calc_local_bh_resistance(flow_rate, temperature) #R_b
        r_v = self.bh_length / (flow_rate * self.fluid.cp(temperature))  # (K/(w/m)) thermal resistance factor


        n = r_v / (r_b * r_a)**(1/2)

        resist_bh_effective_ubt = r_b * n * coth(n)

        return resist_bh_effective_ubt
