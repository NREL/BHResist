from bhr.pipe import Pipe
from math import log

class Borehole:

    def __init__(self, inputs: dict) -> None:

        self.length = inputs["depth"]
        self.pipe = Pipe({**inputs["pipe"], **inputs["fluid"], "length": self.length * 2})
        self.filling = None
        self.shank_space = inputs["shank_space"]
        self.bh_radius = inputs ["borehole-radius"]

        self.theta_1 = self.shank_space / (2 * self.bh_radius)
        self.theta_2 = self.bh_radius / self.pipe.outer_radius
        self.theta_3 = 1 / (2 * self.theta_1 * self.theta_2)

    def calc_bh_average_resistance(self, temperature: float,
                                   flow_rate: float = None,
                                   pipe_resist: float = None) -> float:
        """
        Calculates the average thermal resistance of the borehole using the first-order multipole method.

        Resistance between the fluid in the U-tube(s) to the borehole wall (m-K/W)

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Equation 13

        :param temperature: temperature, Celsius
        :param flow_rate: mass flow rate, kg/s
        :param pipe_resist: pipe thermal resistance, m-K/W
        """

        self.update_beta(temperature, flow_rate, pipe_resist)

        final_term_1 = log(self.theta_2 / (2 * self.theta_1 * (1 - self.theta_1 ** 4) ** self.sigma))

        term_2_num = self.theta_3 ** 2 * (1 - (4 * self.sigma * self.theta_1 ** 4) / (1 - self.theta_1 ** 4)) ** 2
        term_2_den_pt_1 = (1 + self.beta) / (1 - self.beta)
        term_2_den_pt_2 = self.theta_3 ** 2 * (1 + (16 * self.sigma * self.theta_1 ** 4) / (1 - self.theta_1 ** 4) ** 2)
        term_2_den = term_2_den_pt_1 + term_2_den_pt_2
        final_term_2 = term_2_num / term_2_den

        self.resist_bh_ave = (1 / (4 * pi * self.grout.conductivity)) * (self.beta + final_term_1 - final_term_2)

        return self.resist_bh_ave

    def calc_bh_total_internal_resistance(self, temperature: float,
                                          flow_rate: float = None,
                                          pipe_resist: float = None) -> float:
        """
        Calculates the total internal thermal resistance of the borehole using the first-order multipole method.

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Equation 26

        :param temperature: temperature, Celsius
        :param flow_rate: mass flow rate, kg/s
        :param pipe_resist: pipe thermal resistance, m-K/W
        """

        self.update_beta(temperature, flow_rate, pipe_resist)

        term_1_num = (1 + self.theta_1 ** 2) ** self.sigma
        term_1_den = self.theta_3 * (1 - self.theta_1 ** 2) ** self.sigma
        final_term_1 = log(term_1_num / term_1_den)

        term_2_num = self.theta_3 ** 2 * (1 - self.theta_1 ** 4 + 4 * self.sigma * self.theta_1 ** 2) ** 2
        term_2_den_pt_1 = (1 + self.beta) / (1 - self.beta) * (1 - self.theta_1 ** 4) ** 2
        term_2_den_pt_2 = self.theta_3 ** 2 * (1 - self.theta_1 ** 4) ** 2
        term_2_den_pt_3 = 8 * self.sigma * self.theta_1 ** 2 * self.theta_3 ** 2 * (1 + self.theta_1 ** 4)
        term_2_den = term_2_den_pt_1 - term_2_den_pt_2 + term_2_den_pt_3
        final_term_2 = term_2_num / term_2_den

        self.resist_bh_total_internal = 1 / (pi * self.grout.conductivity) * (self.beta + final_term_1 - final_term_2)

        return self.resist_bh_total_internal

    def calc_bh_grout_resistance(self, temperature: float,
                                 flow_rate: float = None,
                                 pipe_resist: float = None) -> float:
        """
        Calculates grout resistance. Use for validation.

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Eq: 3

        :param temperature: temperature, Celsius
        :param flow_rate: mass flow rate, kg/s
        :param pipe_resist: pipe thermal resistance, m-K/W
        """

        self.update_beta(temperature, flow_rate, pipe_resist)

        self.resist_bh_grout = self.calc_bh_average_resistance(temperature, flow_rate,
                                                               pipe_resist) - self.pipe_1.resist_pipe / 2.0
        return self.resist_bh_grout

    def calc_bh_effective_resistance_uhf(self, temperature: float,
                                         flow_rate: float = None,
                                         pipe_resist: float = None) -> float:
        """
        Calculates the effective thermal resistance of the borehole assuming a uniform heat flux.

        Javed, S. & Spitler, J.D. Calculation of Borehole Thermal Resistance. In 'Advances in
        Ground-Source Heat Pump Systems,' pp. 84. Rees, S.J. ed. Cambridge, MA. Elsevier Ltd. 2016.

        Eq: 3-67

        :param temperature: temperature, Celsius
        :param flow_rate: mass flow rate, kg/s
        :param pipe_resist: pipe thermal resistance, m-K/W
        """

        self.update_beta(temperature, flow_rate, pipe_resist)

        self.calc_bh_total_internal_resistance(temperature, flow_rate, pipe_resist)
        self.calc_bh_average_resistance(temperature, flow_rate, pipe_resist)

        pt_1 = 1 / (3 * self.resist_bh_total_internal)
        pt_2 = (self.h / (self.fluid.cp(temperature) * flow_rate)) ** 2
        resist_short_circuiting = pt_1 * pt_2

        self.resist_bh_effective = self.resist_bh_ave + resist_short_circuiting
        return self.resist_bh_effective

    def calc_direct_coupling_resistance(self, temperature: float,
                                        flow_rate: float = None,
                                        pipe_resist: float = None) -> tuple:

        r_a = self.calc_bh_total_internal_resistance(temperature, flow_rate, pipe_resist)
        r_b = self.calc_bh_average_resistance(temperature, flow_rate, pipe_resist)

        r_12 = (4 * r_a * r_b) / (4 * r_b - r_a)

        # reset if negative
        if r_12 < 0:
            r_12 = 70

        self.resist_bh_direct_coupling = r_12
        return self.resist_bh_direct_coupling, r_b

    def update_beta(self, temperature: float, flow_rate: float = None, pipe_resist: float = None) -> float:
        """
        Updates Beta coefficient.

        Javed, S. & Spitler, J.D. Calculation of Borehole Thermal Resistance. In 'Advances in
        Ground-Source Heat Pump Systems,' pp. 84. Rees, S.J. ed. Cambridge, MA. Elsevier Ltd. 2016.

        Eq: 3-47

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Eq: 14

        :param temperature: temperature, Celsius
        :param flow_rate: mass flow rate, kg/s
        :param pipe_resist: pipe thermal resistance, m-K/W
        """

        if flow_rate and pipe_resist:
            # can't set both flow rate and pipe resistance simultaneously
            raise ValueError("'flow_rate' and 'pipe_resist' cannot both be passed.")  # pragma: no cover
        elif flow_rate:
            self.beta = 2 * pi * self.grout.conductivity * self.pipe_1.calc_resist(flow_rate, temperature)
            return self.beta
        elif pipe_resist:
            # setting pipe resistance directly
            # used for validation
            self.pipe_1.resist_pipe = pipe_resist
            self.beta = 2 * pi * self.grout.conductivity * pipe_resist
            return self.beta
        else:
            raise ValueError('Must pass flow rate or a pipe resistance.')  # pragma: no cover
