import unittest

from bhr.fluid import get_fluid


class TestFluid(unittest.TestCase):

    def test_init_ethyl_alcohol(self):
        f = get_fluid(fluid_str="ETHYLALCOHOL", fluid_concentration=0.2)
        self.assertAlmostEqual(f.density(20), 968.9, delta=0.1)

    def test_init_ethylene_glcol(self):
        f = get_fluid(fluid_str="ETHYLENEGLYCOL", fluid_concentration=0.2)
        self.assertAlmostEqual(f.density(20), 1024.1, delta=0.1)

    def test_init_methyl_alcohol(self):
        f = get_fluid(fluid_str="METHYLALCOHOL",  fluid_concentration=0.2)
        self.assertAlmostEqual(f.density(20), 966.7, delta=0.1)

    def test_init_propylene_glcol(self):
        f = get_fluid(fluid_str="PROPYLENEGLYCOL", fluid_concentration=0.2)
        self.assertAlmostEqual(f.density(20), 1014.7, delta=0.1)

    def test_init_water(self):
        f = get_fluid(fluid_str="WATER")
        self.assertAlmostEqual(f.density(20), 998.2, delta=0.1)
