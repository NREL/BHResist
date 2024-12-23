from bhr.pipe import Pipe


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
                 fluid_concentration: float,
                 ):
        self.borehole_diameter = borehole_diameter
        self.grout_conductivity = grout_conductivity
        self.soil_conductivity = soil_conductivity

        self.inner_pipe = Pipe(outer_pipe_outer_diameter, outer_pipe_dimension_ratio, length, outer_pipe_conductivity,
                               fluid_type, fluid_concentration)
        self.outer_pipe = Pipe(inner_pipe_outer_diameter, inner_pipe_dimension_ratio, length, inner_pipe_conductivity,
                               fluid_type, fluid_concentration)

        self.annular_hydraulic_diameter = outer_pipe_outer_diameter - inner_pipe_outer_diameter


    def convective_heat_transfer_coefficients_concentric_annulus(self, flow_rate):


        pass

    def calc_bh_resist(self):
        pass
