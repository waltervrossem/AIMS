#!/usr/bin/env python
# coding: utf-8
# $Id: constants.py
# Author: Daniel R. Reese <daniel.reese@obspm.fr>
# Copyright (C) Daniel R. Reese and contributors
# Copyright license: GNU GPL v3.0
#
#   AIMS is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with AIMS.  If not, see <http://www.gnu.org/licenses/>.
#

r"""
A module which contains the following physical constants:

+------------------------------+--------------------------------------+-------------------------------------+
| **Name of variable**         | **Quantity it describes**            | **Units**                           |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_radius`     | the solar radius                     | :math:`\mathrm{cm}`                 |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_mass`       | the solar mass                       | :math:`\mathrm{g}`                  |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_luminosity` | the solar luminosity                 | :math:`\mathrm{g.cm^2.s^{-3}}`      |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_temperature`| the solar effective temperature      | :math:`\mathrm{K}`                  |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_dnu`        | the solar large frequency separation | :math:`\mathrm{\mu Hz}`             |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_numax`      | the solar frequency at maximum power | :math:`\mathrm{\mu Hz}`             |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_cutoff`     | the solar cutoff frequency           | :math:`\mathrm{\mu Hz}`             |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`G`                | the gravitational constant           | :math:`\mathrm{cm^3.g^{-1}.s^{-2}}` |
+------------------------------+--------------------------------------+-------------------------------------+

.. note::
  These values can be edited according to the latest discoveries.  As
  good practise, it is helpful to include the relevant reference.
"""

__docformat__ = 'restructuredtext'

solar_radius = 6.9598e10
""" the solar radius in :math:`\\mathrm{cm}` """

solar_mass = 1.9892e33
""" the solar mass in :math:`\\mathrm{g}` """

solar_luminosity = 3.8418e33
""" the solar luminosity in :math:`\\mathrm{g.cm^2.s^{-3}}` """

solar_temperature = 5777.0
""" the solar temperature in :math:`\\mathrm{K}` """

solar_dnu = 135.1  # solar delta nu value (in \mu Hz), Huber et al. (2011)
""" the solar large frequency separation in :math:`\\mathrm{\\mu Hz}` """

solar_numax = 3090.0  # solar nu_max value (in \mu Hz), Huber et al. (2011)
""" the solar frequency at maximum power in :math:`\\mathrm{\\mu Hz}` """

solar_cutoff = 5300.0  # Jimenez et al. (2011) (see Balmforth & Gough 1990, Fossat et al. 1992)
""" the solar cut-off frequency separation in :math:`\\mathrm{\\mu Hz}` """

G = 6.67428e-8  # CODATA 2006
""" the gravitational constant in :math:`\\mathrm{cm^3.g^{-1}.s^{-2}}` """
