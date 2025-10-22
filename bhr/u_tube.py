from bhr.pipe import Pipe


class UTube(Pipe):
    def __init__(
        self,
        pipe_outer_diameter: float,
        pipe_dimension_ratio: float,
        length: float,
        shank_space: float,
        pipe_conductivity: float,
        fluid_cp: float = 4182,
        fluid_mu: float = 0.001,
        fluid_rho: float = 998,
        fluid_k: float = 0.598,
    ):
        super().__init__(
            pipe_outer_diameter,
            pipe_dimension_ratio,
            length * 2,
            pipe_conductivity,
            fluid_cp,
            fluid_mu,
            fluid_rho,
            fluid_k,
        )
        self.length = length
        self.shank_space = shank_space
