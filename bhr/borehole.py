from bhr.enums import BoundaryCondition, BoreholeType
from bhr.single_u_borehole import SingleUBorehole
from bhr.double_u_tube import DoubleUTube

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

        if self.bh_type == BoreholeType.SINGLE_U_TUBE:
            bh_diameter = inputs['borehole_diameter']
            pipe_outer_dia = inputs["pipes"]["single-u-tube"]['pipe_outer_diameter']
            dimension_ratio = inputs["pipes"]["single-u-tube"]['pipe_dimension_ratio']
            length = inputs['length']
            shank_space = inputs["pipes"]["single-u-tube"]["shank_space"]
            pipe_conductivity = inputs["pipes"]["single-u-tube"]["pipe_conductivity"]
            grout_conductivity = inputs["grout_conductivity"]
            soil_conductivity = inputs["soil_conductivity"]
            fluid_type = inputs["fluid_type"]

            self.bh = SingleUBorehole(bh_diameter,
                                      pipe_outer_dia,
                                      dimension_ratio,
                                      length,
                                      shank_space,
                                      pipe_conductivity,
                                      grout_conductivity,
                                      soil_conductivity,
                                      fluid_type)
        elif self.bh_type == BoreholeType.DOUBLE_U_TUBE:
            self.bh = DoubleUTube()
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
            return self.bh.calc_bh_effective_resistance_uhf(flow_rate, temperature)

        if self.boundary_condition == BoundaryCondition.UNIFORM_BOREHOLE_WALL_TEMP:
            return self.bh.calc_bh_effective_resistance_ubwt(flow_rate, temperature)
