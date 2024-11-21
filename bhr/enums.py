from enum import Enum, auto


class BoundaryCondition(Enum):
    UNIFORM_BOREHOLE_WALL_TEMP = auto()
    UNIFORM_HEAT_FLUX = auto()


class BoreholeType(Enum):
    SINGLE_U_TUBE = auto()
    DOUBLE_U_TUBE = auto()
    COAXIAL = auto()
