""" This file is intended to calculate the effective borehole
resistance with a double u-pipe configuration. """

from math import pi, log, cosh, sinh

# math functions
ln = log


def coth(x):
    return cosh(x) / sinh(x)


# Input variables needed. To Change in future to real inputs
soil_conductivity = 3  # thermal conductivity of the ground W/(m-K) #
grout_conductivity = 1.5  # thermal conductivity of the borehole/grout W/(m-K) #
rb = 57.5  # radius of borehole (mm) #
rp = 16  # radius of pipe (mm) #
length = 200  # depth of borehole (m) #
Rp = 0.05  # total fluid-to-pipe thermal resistance for a single pipe K/(W/m)
Vf = 0.75 / 3600  # (m^3/s) volumetric flow rate for one pipe, one direction
cf = 4180  # specific heat of circulating fluid (J/kg)/K
pf = 997  # (kg/m3) density of circulating fluid
rc = 22.63  # (mm) radial distance between centers of symmetrically placed pipes and borehole center

print("rc = %.2f mm" % rc)

# check to see if rc is rp*sqrt(2) < rc < rb-rp
# (add code here)

# -local derived variables-
Rv = length / (pf * cf * Vf)  # (K/(w/m)) thermal resistance factor
Ppc = rp ** 2 / (4 * rc ** 2)  # dimensionless parameter
Pc = rc ** 2 / (rb ** 8 - rc ** 8) ** (1 / 4)  # dimensionless parameter
Pb = rb ** 2 / (rb ** 8 - rc ** 8) ** (1 / 4)  # dimensionless parameter
beta = 2 * pi * grout_conductivity * Rp  # dimensionless thermal resistance of 1 u-pipe leg
b1 = (1 - beta) / (1 + beta)  # dimensionless parameter
sigma = (grout_conductivity - soil_conductivity) / (
            grout_conductivity + soil_conductivity)  # thermal conductivity ratio, dimensionless

print("""
 Rv = %.2E
 Ppc = %.2E
 Pc = %.2E
 Pb = %.2E
 beta = %.2E
 b1 = %.2E
 sigma = %.2E \n""" % (Rv, Ppc, Pc, Pb, beta, b1, sigma))

# --Borehole resistance, 0th order [K/(W/m)]--
Rb0 = Rp / 4 + 1 / (8 * pi * grout_conductivity) * ((ln(rb ** 4 / (4 * rp * rc ** 3))) +
                                                    sigma * ln(rb ** 8 / (rb ** 8 - rc ** 8)))

# --Borehole resistance, 1st order [K/(W/m)]--
Rb1 = Rb0 - 1 / (8 * pi * grout_conductivity) * (b1 * Ppc * (3 - 8 * sigma * Pc ** 4) ** 2
                                                 ) / (1 + b1 * Ppc * (5 + 64 * sigma * Pc ** 4 * Pb ** 4))

print("Rb0 = %.3E" % Rb0)
print("Rb1 = %.3E" % Rb1)

# --Internal thermal resistance for diagonal inlet pipes [K/(W/m)]--

# 0th order
Rad0 = 2 * Rp + 2 / (2 * pi * grout_conductivity) * (ln(rc / rp) +
                                                     sigma * ln((rb ** 4 + rc ** 4) / (rb ** 4 - rc ** 4)))

# 1st order
Rad1 = Rad0 - 2 / (2 * pi * grout_conductivity) * (b1 * Ppc * (1 + 8 * sigma * Pc ** 2 * Pb ** 2) ** 2
                                                   ) / (
                   1 - b1 * Ppc * (3 - 32 * sigma * (Pc ** 2 * Pb ** 6 + Pc ** 6 * Pb ** 2)))

print("Rad0 = %.3E" % Rad0)
print("Rad1 = %.3E" % Rad1)

# --Internal thermal resistance for adjacent inlet pipes [K/(W/m)]--
# 0th order
Raa0 = 2 * Rp + 2 / (2 * pi * grout_conductivity) * (ln(2 * rc / rp) +
                                                     sigma * ln((rb ** 2 + rc ** 2) / (rb ** 2 - rc ** 2)))

# 1st order
M11 = 1 + 16 * b1 * sigma * Ppc * (3 * Pc ** 3 * Pb ** 5 + Pc ** 7 * Pb)  # matrix variable
M22 = -1 - 16 * b1 * sigma * Ppc * (Pc * Pb ** 7 + 3 * Pc ** 5 * Pb ** 3)  # matrix variable
M21 = b1 * Ppc  # matrix variable
V1 = 1 - 8 * sigma * Pc ** 3 * Pb  # vector variable
V2 = 3 + 8 * sigma * Pc * Pb ** 3  # vector variable

Raa1 = Raa0 + 2 / (2 * pi * grout_conductivity) * b1 * Ppc / 2 * (V2 ** 2 * M11 - 2 * V1 * V2 * M21 -
                                                                  V1 ** 2 * M22) / (M11 * M22 + M21 ** 2)

print("Raa0 = %.3E" % Raa0)
print("Raa1 = %.3E\n" % Raa1)

# ---Uniform heat flux effective borehole resistance, 1st order---

# diagonal
Rb_eff_d_UHF = Rb1 + Rv ** 2 / (6 * Rad1)
# adjacent
Rb_eff_a_UHF = Rb1 + Rv ** 2 / (6 * Raa1)

print("Rb_eff_d_UHF = %.3E" % Rb_eff_d_UHF)
print("Ra_eff_a_UHF = %.3E" % Rb_eff_a_UHF)

# ---Uniform wall temperature, 1st order---

# diagonal
n_d = Rv / (2 * Rb1 * Rad1) ** (1 / 2)
Rb_eff_d_UWT = Rb1 * n_d * coth(n_d)

# adjacent
n_a = Rv / (2 * Rb1 * Raa1) ** (1 / 2)
Rb_eff_a_UWT = Rb1 * n_a * coth(n_a)

print("Rb_eff_d_UWT = %.3E" % Rb_eff_d_UWT)
print("Ra_eff_a_UWT = %.3E\n" % Rb_eff_a_UWT)

# ----Final effective borehole resistance [K/(W/m)]----
# diagonal, average of UWT and UHF
Rb_eff_d = (Rb_eff_d_UHF + Rb_eff_d_UWT) / 2
# adjacent, average of UWT and UHF
Rb_eff_a = (Rb_eff_a_UHF + Rb_eff_a_UWT) / 2

print("Rb_eff_d = %.3E [K/(W/m)]" % Rb_eff_d)
print("Ra_eff_a = %.3E [K/(W/m)]" % Rb_eff_a)
