import unittest

from bhr.fluid import get_fluid


class TestFluid(unittest.TestCase):
    def test_init_ethyl_alcohol(self):
        f = get_fluid(fluid_type="ETHYLALCOHOL", fluid_concentration=0.2)
        self.assertAlmostEqual(f.density(20), 968.9, delta=0.1)

    def test_init_ethylene_glycol(self):
        f = get_fluid(fluid_type="ETHYLENEGLYCOL", fluid_concentration=0.2)
        self.assertAlmostEqual(f.density(20), 1024.1, delta=0.1)

    def test_init_methyl_alcohol(self):
        f = get_fluid(fluid_type="METHYLALCOHOL", fluid_concentration=0.2)
        self.assertAlmostEqual(f.density(20), 966.7, delta=0.1)

    def test_init_propylene_glycol(self):
        f = get_fluid(fluid_type="PROPYLENEGLYCOL", fluid_concentration=0.2)
        self.assertAlmostEqual(f.density(20), 1014.7, delta=0.1)

    def test_init_water(self):
        f = get_fluid(fluid_type="WATER")
        self.assertAlmostEqual(f.density(20), 998.2, delta=0.1)
