from bhr.double_u_borehole import DoubleUTube
from bhr.enums import BoundaryCondition, BoreholeType
from bhr.single_u_borehole import SingleUBorehole


class Borehole:

    def __init__(self):
        self.bh_type = None
        self.boundary_condition = None
        self.bh = None

    def init_from_dict(self, inputs: dict):

        bh_type_str = inputs["borehole_type"].upper()
        if bh_type_str == BoreholeType.SINGLE_U_TUBE.name:
            self.bh_type = BoreholeType.SINGLE_U_TUBE
        elif bh_type_str == BoreholeType.DOUBLE_U_TUBE.name:
            self.bh_type = BoreholeType.DOUBLE_U_TUBE
        elif bh_type_str == BoreholeType.COAXIAL.name:
            self.bh_type = BoreholeType.COAXIAL
        else:
            raise LookupError(f"borehole_type \"{bh_type_str}\" not supported")

        bc_str = inputs["boundary_condition"].upper()
        if bc_str == BoundaryCondition.UNIFORM_HEAT_FLUX.name:
            self.boundary_condition = BoundaryCondition.UNIFORM_HEAT_FLUX
        elif bc_str == BoundaryCondition.UNIFORM_BOREHOLE_WALL_TEMP.name:
            self.boundary_condition = BoundaryCondition.UNIFORM_BOREHOLE_WALL_TEMP
        else:
            raise LookupError(f"boundary_condition \"{bc_str}\" not supported")

        bh_diameter = inputs['borehole_diameter']
        length = inputs['length']
        grout_conductivity = inputs["grout_conductivity"]
        soil_conductivity = inputs["soil_conductivity"]
        fluid_type = inputs["fluid_type"]
        fluid_concentration = inputs["fluid_concentration"]

        if self.bh_type == BoreholeType.SINGLE_U_TUBE:
            pipe_outer_dia_single = inputs["single_u_tube"]['pipe_outer_diameter']
            dimension_ratio_single = inputs["single_u_tube"]['pipe_dimension_ratio']
            shank_space_single = inputs["single_u_tube"]["shank_space"]
            pipe_conductivity_single = inputs["single_u_tube"]["pipe_conductivity"]

            self.bh = SingleUBorehole(bh_diameter,
                                      pipe_outer_dia_single,
                                      dimension_ratio_single,
                                      length,
                                      shank_space_single,
                                      pipe_conductivity_single,
                                      grout_conductivity,
                                      soil_conductivity,
                                      fluid_type,
                                      fluid_concentration)

        elif self.bh_type == BoreholeType.DOUBLE_U_TUBE:

            pipe_outer_dia_double = inputs["double_u_tube"]['pipe_outer_diameter']
            dimension_ratio_double = inputs["double_u_tube"]['pipe_dimension_ratio']
            shank_space_double = inputs["double_u_tube"]["shank_space"]
            pipe_conductivity_double = inputs["double_u_tube"]["pipe_conductivity"]
            pipe_inlet_arrangement = inputs["double_u_tube"]["pipe_inlet_arrangement"]

            self.bh = DoubleUTube(bh_diameter,
                                  pipe_outer_dia_double,
                                  dimension_ratio_double,
                                  length,
                                  shank_space_double,
                                  pipe_conductivity_double,
                                  pipe_inlet_arrangement,
                                  grout_conductivity,
                                  soil_conductivity,
                                  fluid_type,
                                  fluid_concentration)

        elif self.bh_type == BoreholeType.COAXIAL:
            pass
        else:
            raise NotImplementedError(f"bh_type \"{self.bh_type.name}\" not implemented")

    def init_from_file(self, file_path):
        # setup input dict
        # call init_from_dict()
        pass

    def init_from_params(self):
        pass

    def _setup(self):
        pass

    def calc_bh_resist(self, flow_rate, temperature):

        if self.bh is None:
            pass
            # assert error and message...

        if self.boundary_condition == BoundaryCondition.UNIFORM_HEAT_FLUX:
            return self.bh.calc_effective_bh_resistance_uhf(flow_rate, temperature)

        if self.boundary_condition == BoundaryCondition.UNIFORM_BOREHOLE_WALL_TEMP:
            return self.bh.calc_effective_bh_resistance_ubwt(flow_rate, temperature)
