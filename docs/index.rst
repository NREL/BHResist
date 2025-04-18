.. BH Resist documentation master file, created by
   sphinx-quickstart on Wed Feb 19 11:27:40 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

BHResist Documentation
=======================

A pure python library for computing thermal resistance within single-u, double-u,
and coaxial grouted borehole heat exchangers. For single and double u-tube configurations,
the methods use the 1st-order closed-form multipole approximations, which typically
produces results with less than 1% error when compared to the 10th-order multipole method.
Coaxial borehole methods apply a simple 1D resistance network method.

This is intended to be a lightweight library that can be easily imported into any other Python tool,
with no bulky dependencies.

References
----------

Hellström, G. 1991. "Ground Heat Storage: Thermal Analyses of Duct Storage Systems." PhD dissertation.
Department of Mathematical Physics, University of Lund, Sweden.

Grundmann, R.M. 2016 "Improved design methods for ground heat exchangers." Master’s thesis, Oklahoma State University.

Javed, S. and J.D. Spitler. 2016. "Calculation of borehole thermal resistance." In *Advances in Ground-Source
Heat Pump Systems*. Ed. S.J. Rees. Woodhead Publishing. https://doi.org/10.1016/B978-0-08-100311-4.00003-0

Javed, S., and J.D. Spitler. 2017. "Accuracy of borehole thermal resistance calculation methods for grouted
single u-tube ground heat exchangers." *Applied Energy,* 187:790-806. https://doi.org/10.1016/j.apenergy.2016.11.079

Claesson, J., and S. Javed. 2019. "Explicit multipole formulas and thermal network models for calculating
thermal resistances of double U-pipe borehole heat exchangers." *Science and Technology for the Built
Environment,* 25(8) pp. 980–992. https://doi.org/10.1080/23744731.2019.1620565

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   prog_usage
   borehole
   single_bh
   double_bh
   coaxial_bh
