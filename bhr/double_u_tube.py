from bhr.u_tube import UTube

from math import pi

class DoubleUTube(UTube):

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

        # static parameters
        self.borehole_diameter = borehole_diameter
        self.grout_conductivity = grout_conductivity
        self.soil_conductivity = soil_conductivity
        self.sigma = (self.grout_conductivity - self.soil_conductivity) / (
                self.grout_conductivity + self.soil_conductivity)

        # non-static parameters
        self.beta = None

    def update_beta(self,
                    flow_rate: float = None,
                    temperature: float = None,
                    pipe_resist=None):
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
        :param pipe_resist: pipe conduction and convection resistance, K/(W/m). only used for testing
        """

        if pipe_resist:
            self.pipe_resist = pipe_resist
        else:
            self.pipe_resist = self.calc_pipe_resist(flow_rate, temperature)

        self.beta = 2 * pi * self.grout_conductivity * self.pipe_resist

    def calc_bh_resist(self, flow_rate, temperature):
        pass

    def calc_internal_bh_resist_pipe(self, flow_rate, temperature, pipe_config):

        if pipe_config == "DIAGONAL":
            pass
        elif pipe_config == "ADJACENT":
            pass
        else:
            assert False

    def calc_effective_bh_resist_uhf(self, flow_rate, temperature):
        return 0.2

    def calc_effective_bh_resist_ubwt(self, flow_rate, temperature):
        return 0.2

    def calc_effective_bh_resist_ave(self, flow_rate, temperature):
        rb_uhf = self.calc_effective_bh_resist_uhf(flow_rate, temperature)
        rb_ubwt = self.calc_effective_bh_resist_ubwt(flow_rate, temperature)
        return (rb_uhf + rb_ubwt) / 2.0
