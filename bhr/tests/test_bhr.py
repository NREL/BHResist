# import unittest
#
# from bhr.single_u_borehole import SingleUBorehole
#
#
# class TestBHR(unittest.TestCase):
#
#     def test_init_single_u_tube_bh(self):
#         inputs = {
#             "fluid": {
#                 "fluid_type": "PROPYLENEGLYCOL",
#                 "fluid_concentration": 0.2,
#             },
#             "pipes": {
#                 "pipe-type": "u-tube",
#                 "u-tube": {
#                     "outer_diameter": 0.042,
#                     "num-u_tubes": 1,
#                     "dimension-ratio": 11,
#                     "conductivity": 0.4,
#                     "shank_space": 0.01,
#                     "hydraulic_configuration": "PARALLEL"  # or SERIES
#                 },
#                 "coaxial": {
#                     "inner_pipe": {
#                         "outer_diameter": 0.042,
#                         "dimension_ratio": 11,
#                         "conductivity": 0.4
#                     },
#                     "outer-pipe": {
#                         "outer_diameter": 0.096,
#                         "dimension_ratio": 11,
#                         "conductivity": 0.4
#                     }
#                 }
#             },
#             "filling": {
#                 "type": "GROUT",
#                 "conductivity": 1.2,
#                 "heat_capacity": 1600000
#             },
#             "depth": 100,
#             "borehole_radius": 0.14
#         }
#
#         bh = SingleUBorehole(inputs)
#
#         # only pass flow rate, so pipe resistance should be computed in the process of this call
#         bh.calc_bh_effective_resistance_uhf(temperature=20, flow_rate=0.5)
