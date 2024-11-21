from bhr.u_tube import UTube


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

        # non-static parameters

    def calc_effective_bh_resist_uhf(self, flow_rate, temperature):
        return 0.2

    def calc_effective_bh_resist_ubwt(self, flow_rate, temperature):
        return 0.2
