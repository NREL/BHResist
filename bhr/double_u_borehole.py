from math import pi, log, cosh, sinh

from bhr.u_tube import UTube

ln = log


def coth(x):
    return cosh(x) / sinh(x)


class DoubleUTube(UTube):

    def __init__(self,
                 borehole_diameter: float,
                 pipe_outer_diameter: float,
                 pipe_dimension_ratio: float,
                 length: float,
                 shank_space: float,
                 #  pipe_resist: float,
                 pipe_conductivity: float,
                 pipe_config: str,
                 grout_conductivity: float,
                 soil_conductivity: float,
                 fluid_type: str,
                 fluid_concentration: float = 0):

        super().__init__(pipe_outer_diameter, pipe_dimension_ratio, length, shank_space, pipe_conductivity,
                         fluid_type, fluid_concentration)

        # static parameters
        self.grout_conductivity = grout_conductivity
        self.borehole_radius = borehole_diameter / 2  # radius of borehole (mm) rb
        self.pipe_radius = pipe_outer_diameter / 2  # pipe outer radius (mm) rp
        self.pipe_config = pipe_config
        self.length = length  # length of borehole (m) H
        self.grout_conductivity = grout_conductivity
        self.soil_conductivity = soil_conductivity
        # self.pipe_resist = pipe_resist #total fluid-to-pipe resistance for a single pipe K/(W/m) Rp
        # self.fluid_density = fluid_density #(kg/m3) density of circulating fluid pf
        # self.Vf = mass_flow_rate / fluid_density #(m^3/s) volumetric flow rate for one pipe, one direction Vf
        # self.cf = fluid_specific_heat # specific heat of circulating fluid (J/kg)/K cf
        self.shank_space = shank_space  # (mm) radial distance between centers of symmetrically placed pipes and borehole center rc
        self.sigma = (self.grout_conductivity - self.soil_conductivity) / (
                self.grout_conductivity + self.soil_conductivity)

        # -derived variables-
        self.Ppc = self.pipe_radius ** 2 / (4 * self.shank_space ** 2)  # dimensionless parameter
        self.Pc = self.shank_space ** 2 / (self.borehole_radius ** 8 - self.shank_space ** 8) ** (
                    1 / 4)  # dimensionless parameter
        self.Pb = self.borehole_radius ** 2 / (self.borehole_radius ** 8 - self.shank_space ** 8) ** (
                    1 / 4)  # dimensionless parameter

        self.sigma = (grout_conductivity - soil_conductivity) / (
                    grout_conductivity + soil_conductivity)  # thermal conductivity ratio, dimensionless

        # non-static parameters
        self.beta = None
        self.Rb1 = None
        self.Ra1 = None
        self.Rb_eff_UHF = None
        self.Rb_eff_UWT = None
        self.eff_bhr_ave = None

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
        self.b1 = (1 - self.beta) / (1 + self.beta)  # dimensionless parameter

    def calc_bh_resist(self, flow_rate, temperature):
        """
        Calculates tube-to-borehole resistance (aka local borehole resistnance) .

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
         the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

         eqns.13 & 14
         """

        self.update_beta(flow_rate, temperature)

        # --Borehole resistance, 0th order [K/(W/m)]--
        Rb0 = self.pipe_resist / 4 + 1 / (8 * pi * self.grout_conductivity) * (
                (ln(self.borehole_radius ** 4 / (4 * self.pipe_radius * self.shank_space ** 3))) +
                self.sigma * ln(self.borehole_radius ** 8 / (self.borehole_radius ** 8 - self.shank_space ** 8)))

        # --Borehole resistance, 1st order [K/(W/m)]--
        self.Rb1 = Rb0 - 1 / (8 * pi * self.grout_conductivity) * (
                    self.b1 * self.Ppc * (3 - 8 * self.sigma * self.Pc ** 4) ** 2
                    ) / (1 + self.b1 * self.Ppc * (5 + 64 * self.sigma * self.Pc ** 4 * self.Pb ** 4))

        #print("Rb0 = %.3E" % Rb0)
        #print("Rb1 = %.3E" % self.Rb1)

        return self.Rb1

    def calc_internal_bh_resist_pipe(self, flow_rate, temperature):
        """
        Calculates tube-to-tube resistance (aka internal resistance).

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
         the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

         eqns. 18, 19, 22, 23
         """

        self.update_beta(flow_rate, temperature)

        if self.pipe_config == "DIAGONAL":
            # 0th order
            Ra0 = 2 * self.pipe_resist + 2 / (2 * pi * self.grout_conductivity) * (
                        ln(self.shank_space / self.pipe_radius) +
                        self.sigma * ln((self.borehole_radius ** 4 + self.shank_space ** 4) / (
                            self.borehole_radius ** 4 - self.shank_space ** 4)))

            # 1st order
            self.Ra1 = Ra0 - 2 / (2 * pi * self.grout_conductivity) * (
                        self.b1 * self.Ppc * (1 + 8 * self.sigma * self.Pc ** 2 * self.Pb ** 2) ** 2
                        ) / (1 - self.b1 * self.Ppc * (
                        3 - 32 * self.sigma * (self.Pc ** 2 * self.Pb ** 6 + self.Pc ** 6 * self.Pb ** 2)))

            #print("Rad0 = %.3E" % Ra0)
            #print("Rad1 = %.3E" % self.Ra1)

            return self.Ra1

        elif self.pipe_config == "ADJACENT":
            # 0th order
            Ra0 = 2 * self.pipe_resist + 2 / (2 * pi * self.grout_conductivity) * (
                        ln(2 * self.shank_space / self.pipe_radius) +
                        self.sigma * ln((self.borehole_radius ** 2 + self.shank_space ** 2) / (
                            self.borehole_radius ** 2 - self.shank_space ** 2)))

            # 1st order
            M11 = 1 + 16 * self.b1 * self.sigma * self.Ppc * (
                        3 * self.Pc ** 3 * self.Pb ** 5 + self.Pc ** 7 * self.Pb)  # matrix variable
            M22 = -1 - 16 * self.b1 * self.sigma * self.Ppc * (
                        self.Pc * self.Pb ** 7 + 3 * self.Pc ** 5 * self.Pb ** 3)  # matrix variable
            M21 = self.b1 * self.Ppc  # matrix variable
            V1 = 1 - 8 * self.sigma * self.Pc ** 3 * self.Pb  # vector variable
            V2 = 3 + 8 * self.sigma * self.Pc * self.Pb ** 3  # vector variable

            self.Ra1 = Ra0 + 2 / (2 * pi * self.grout_conductivity) * self.b1 * self.Ppc / 2 * (
                        V2 ** 2 * M11 - 2 * V1 * V2 * M21 -
                        V1 ** 2 * M22) / (M11 * M22 + M21 ** 2)

            return self.Ra1
        else:
            assert False

    def calc_effective_bh_resist_uhf(self, flow_rate, temperature):
        """
       Calculates effective borehole resistance for uniform heat flux along the borehole.

       Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
       and Thermal Network Models for Calculating Thermal Resistances of
       Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

        eq. 44
        """
        Rv = self.length / ( self.fluid.cp(temperature) * flow_rate)  # (K/(w/m)) thermal resistance factor

        self.Rb_eff_UHF = self.Rb1 + Rv ** 2 / (6 * self.Ra1)
        # print("Rb_eff_d_UHF = %.3E" % self.Rb_eff_UHF)

        return self.Rb_eff_UHF

    def calc_effective_bh_resist_ubwt(self, flow_rate, temperature):
        """
        Calculates effective borehole resistance for uniform borehole wall temperature.

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
         the Built Environment 25 (8): 980–92. doi:10.1080/23744731.2019.1620565.

         eq. 46
         """
        Rv = self.length / (flow_rate * self.fluid.cp(temperature))  # (K/(w/m)) thermal resistance factor
        n = Rv / (2 * self.Rb1 * self.Ra1) ** (1 / 2)
        self.Rb_eff_UWT = self.Rb1 * n * coth(n)

        #print("Ra_eff_UWT = %.3E\n" % self.Rb_eff_UWT)

        return self.Rb_eff_UWT

    def calc_effective_bh_resist_ave(self):
        """ Averages effective borehole resistance results for the two boundry conditions
        Uniform Heat Flux and Uniform Borehole Wall temperature
        """
        self.eff_bhr_ave = (self.Rb_eff_UHF + self.Rb_eff_UWT) / 2.0
        return self.eff_bhr_ave
