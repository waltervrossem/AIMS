#!/usr/bin/env python
# coding: utf-8
# $Id: model.py
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

"""
A module which contains various classes relevant to the grid of models:

- :py:class:`Model`: a model
- :py:class:`Track`: an evolutionary track
- :py:class:`Model_grid`: a grid of models

These different classes allow the program to store a grid of models and perform
a number of operations, such as:

- retrieving model properties
- interpolate within the grid models
- sort the models within a given evolutionary track
- ...

"""

__docformat__ = 'restructuredtext'

import os
import sys

# AIMS configuration:
if 'AIMS_configure.py' in os.listdir(os.getcwd()):
    sys.path = [os.getcwd(), *sys.path]
import AIMS_configure as config

print(f'AIMS_configure.py path in model.py: {config.__file__}')
# packages from within AIMS:
import constants
import utilities
import aims_fortran

# various packages needed for AIMS
import math
import numpy as np
import random
import matplotlib

if (config.backend is not None):
    matplotlib.use(config.backend)
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from bisect import bisect_right
import h5py

# maximum filename length
# filename_dtype = "U1000"
filename_dtype = "str"

# Strings for age parameters
age_str = "Age"
""" Name of the physical age parameter used in the interpolation """

age_adim_str = "Age_adim"
""" Name of the non-dimensional age parameter used in the interpolation """

# Global parameters:
tol = 1e-6
""" tolerance level for points slightly outside the grid """

eps = 1e-6
""" relative tolerance on parameters used for setting up evolutionary tracks """

log0 = -1e150
""" value which is returned when log(0) is calculated (rather than causing an error) """

ntype = np.int16
""" type used for the n values """

ltype = np.int8  # int8 ranges from -128 to 127
""" type used for the l values """

ftype = np.float64
"""  type used for the frequencies """

gtype = np.float64
""" type used for grid data """

modetype = [('n', ntype), ('l', ltype), ('freq', ftype), ('inertia', ftype)]
""" structure for modes """

# quantities related to user-defined parameters:
user_params_index = {}
""" dictionary which will supply the appropriate index for the user-defined parameters"""

user_params_latex = {}
""" dictionary which will supply the appropriate latex name for the user-defined parameters"""

# integer indices for various global quantities which will be stored in a np array:
nglb = 9 + len(config.user_params)
""" total number of global quantities in a model (see :py:data:`Model.glb`)."""

nlin = 6 + len(config.user_params)
"""
total number of global quantities which are interpolated in a linear way (see
:py:func:`combine_models`).  These quantities are numbered 0:nlin-1
"""

iage_adim = 0
""" index of the parameter corresponding to dimensionless age in the :py:data:`Model.glb` array """

iage = 1
""" index of the parameter corresponding to age in the :py:data:`Model.glb` array """

imass = 2
""" index of the parameter corresponding to mass in the :py:data:`Model.glb` array """

itemperature = 3
""" index of the parameter corresponding to temperature in the :py:data:`Model.glb` array """

iz0 = 4
"""
index of the parameter corresponding to the initial metallicity
the :py:data:`Model.glb` array
"""

ix0 = 5
"""
index of the parameter corresponding to the initial hydrogen content
in the :py:data:`Model.glb` array
"""

ifreq_ref = 6 + len(config.user_params)
"""
index of the parameter corresponding to the reference frequency
(used to non-dimensionalise the pulsation frequencies of the model)
in the :py:data:`Model.glb` array
"""

iradius = 7 + len(config.user_params)
""" index of the parameter corresponding to radius in the :py:data:`Model.glb` array """

iluminosity = 8 + len(config.user_params)
""" index of the parameter corresponding to luminosity in the :py:data:`Model.glb` array """


def string_to_latex(string, prefix="", postfix=""):
    """
    Return a fancy latex name for an input string.
    
    :param string: string that indicates for which parameter we're seeking a latex name
    :param prefix: optional prefix to add to the string
    :param postfix: optional postfix to add to the string

    :type string: string
    :type prefix: string
    :type postfix: string

    :return: a fancy latex name
    :rtype: string

    .. note::
      This also works for the names of the amplitude parameters for surface corrections.
    """

    if (string.startswith("log_")):
        return string_to_latex(string[4:], prefix + r"\log_{10}\left(", r"\right)" + postfix)
    if (string.startswith("ln_")):
        return string_to_latex(string[3:], prefix + r"\ln\left(", r"\right)" + postfix)
    if (string.startswith("exp_")):
        return string_to_latex(string[4:], prefix + r"\exp\left(", r"\right)" + postfix)
    if (string == "Mass"):
        return r'Mass, $%sM/M_{\odot}%s$' % (prefix, postfix)
    if (string == "Radius"):
        return r'Radius, $%sR/R_{\odot}%s$' % (prefix, postfix)
    if (string == "Luminosity"):
        return r'Luminosity, $%sL/L_{\odot}%s$' % (prefix, postfix)
    if (string == "Z"):
        return r'Metallicity, $%sZ_0%s$' % (prefix, postfix)
    if (string == "Y"):
        return r'Helium content, $%sY_0%s$' % (prefix, postfix)
    if (string == "X"):
        return r'Hydrogen content, $%sX_0%s$' % (prefix, postfix)
    if (string == "Ys"):
        return r'Helium content, $%sY_s%s$' % (prefix, postfix)
    if (string == "Yc"):
        return r'Central helium, $%sY_s%s$' % (prefix, postfix)
    if (string == "zsx_s"):
        return r'Ratio, $%s(Z/X)_s%s$' % (prefix, postfix)
    if (string == "zsx_0"):
        return r'Ratio, $%s(Z/X)_0%s$' % (prefix, postfix)
    if (string == "Fe_H"):
        return r'Iron content, $%s\mathrm{[Fe/H]}%s$' % (prefix, postfix)
    if (string == "Fe_H0"):
        return r'Initial iron content, $%s\mathrm{[Fe/H]}%s$' % (prefix, postfix)
    if (string == "M_H"):
        return r'Metal content, $%s\mathrm{[M/H]}%s$' % (prefix, postfix)
    if (string == "M_H0"):
        return r'Initial metal content, $%s\mathrm{[M/H]}_0%s$' % (prefix, postfix)
    if (string == "dY_dZ"):
        return r'Enrichment law, $%s\mathrm{dY/dZ}%s$' % (prefix, postfix)
    if (string == age_str):
        return r'Age (in Myrs), $%st%s$' % (prefix, postfix)
    if (string == age_adim_str):
        return r'Age parameter, $%s\tau%s$' % (prefix, postfix)
    if (string == "Teff"):
        return r'$%sT_{\mathrm{eff}}%s$ (in K)' % (prefix, postfix)
    if (string == "Dnu"):
        return r'Large separation (in $\mu$Hz), $%s\Delta \nu%s$' % (prefix, postfix)
    if (string == "Dnu_c"):
        return r'Corrected Large separation (in $\mu$Hz), $%s\Delta \nu%s$' % (prefix, postfix)
    if (string == "numax"):
        return r'$%s\nu_{\mathrm{max}}%s$ (in $\mu$Hz)' % (prefix, postfix)
    if (string == "Rho"):
        return r'Density (in g.cm$^{-3}$), $%s\rho%s$' % (prefix, postfix)
    if (string == "g"):
        return r'Surface gravity (in cm.s$^{-2}$), $%sg%s$' % (prefix, postfix)
    if (string == "Theta_BCZ"):
        return r'BCZ normalised acoustic depth, $%s\tau_{\mathrm{BCZ}}/2\tau%s$' % (
            prefix, postfix)
    if (string == "Theta_He"):
        return r'HeII normalised acoustic depth, $%s\tau_{\mathrm{He}}/2\tau%s$' % (
            prefix, postfix)
    if (string == "A_surf"):
        return r'$%sA_{\mathrm{surf}}%s$' % (prefix, postfix)
    if (string == "A3_surf"):
        return r'$%sA_3^{\mathrm{surf}}%s$' % (prefix, postfix)
    if (string == "Am1_surf"):
        return r'$%sA_{-1}^{\mathrm{surf}}%s$' % (prefix, postfix)
    if (string == "alpha_surf"):
        return r'$%s\alpha_{\mathrm{surf}}%s$' % (prefix, postfix)
    if (string == "b_Kjeldsen2008"):
        return r'$%sb_{\mathrm{Kjeldsen\,et\,al.\,(2008)}}%s$' % (prefix, postfix)
    if (string == "beta_Sonoi2015"):
        return r'$%s\beta_{\mathrm{Sonoi\,et\,al.\,(2015)}}%s$' % (prefix, postfix)

    # raises a KeyError if `string` isn't a valid key
    return user_params_latex[string] % (prefix, postfix)


def get_surface_parameter_names(surface_option):
    """
    Return the relevant parameter names for a given surface correction option.

    :param surface_option: specifies the type of surface correction.
    :type surface_option: string

    :return: names for the surface parameters
    :rtype: tuple of strings
    """

    if (surface_option is None):
        return ()
    if (surface_option == "Kjeldsen2008"):
        return ("A_surf",)
    if (surface_option == "Kjeldsen2008_scaling"):
        return ("A_surf",)
    if (surface_option == "Kjeldsen2008_2"):
        return ("A_surf", "b_Kjeldsen2008")
    if (surface_option == "Ball2014"):
        return ("A3_surf",)
    if (surface_option == "Ball2014_2"):
        return ("A3_surf", "Am1_surf")
    if (surface_option == "Sonoi2015"):
        return ("alpha_surf",)
    if (surface_option == "Sonoi2015_scaling"):
        return ("alpha_surf",)
    if (surface_option == "Sonoi2015_2"):
        return ("alpha_surf", "beta_Sonoi2015")
    raise ValueError("Unknown surface correction: " + surface_option)

def get_aFe(track):
    """Get [alpha/Fe] from a track."""
    if isinstance(config.alpha_Fe_param, str):
        if config.alpha_Fe_param in track.grid_params:
            i_aFe = track.grid_params.index(config.alpha_Fe_param)
            aFe = track.params[i_aFe]
        else:
            aFe = None
    else:
        aFe = config.alpha_Fe_param
    return aFe

class Model:
    """A class which contains a stellar model, including classical and seismic information."""

    def __init__(self, _glb, _name=None, _modes=None, aFe=None):
        """
        :param _glb: 1D array of global parameters for this model.  Its dimension should
          be greater or equal to :py:data:`nglb`
        :param _name: name of the model (typically the second part of its path)
        :param _modes: list of modes in the form of tuples
          (n,l,freq,inertia) which will be appended to the set of modes in the model.
        :param aFe: [alpha/Fe] used to calculate A_FeH correction factor. Can be a float
          or string or None. If a string, it must be a grid parameter.

        :type _glb: np.array
        :type _name: string
        :type _modes: list of (int, int, float, float)
        """

        # check inputs
        assert (_glb[imass] >= 0.0), "A star cannot have a negative mass!"
        assert (_glb[iradius] >= 0.0), "A star cannot have a negative radius!"
        assert (_glb[iluminosity] >= 0.0), "A star cannot have a negative luminosity!"
        assert (_glb[iz0] >= 0.0), "A star cannot have a negative metallicity!"
        assert (_glb[ix0] >= 0.0), "A star cannot have a negative hydrogen abundance!"
        assert (_glb[iage] >= 0.0), "A star cannot have a negative age!"
        assert (_glb[itemperature] >= 0.0), "A star cannot have a negative temperature!"

        self.name = _name
        """Name of the model, typically the second part of its path"""

        self.glb = _glb
        """Array which will contain various global quantities"""

        self.glb[ifreq_ref] = 5e5 * math.sqrt(constants.G * self.glb[imass] / self.glb[iradius] ** 3) / math.pi
        """Characteristic frequency of the model in :math:`\\mathrm{cyclic \\, \\mu Hz}`"""

        """array containing the modes (n, l, freq, inertia)"""
        if _modes is not None:
            self.modes = np.array(_modes, dtype=modetype)
        else:
            self.modes = np.empty([0], dtype=modetype)

        if config.alpha_Fe_param is None or aFe is None:
            self.A_FeH = 1.0
        else:
            # These coefficients are calculated using AIMS/etc/calc_A_FeH_coeffs.py
            if config.use_Asplund_A_FeH:
                # Using Asplund 2009 solar mixture
                a = 0.6386992237569251
                b = 0.3613007762430738
            else:
                # Using gs98 solar mixture
                a = 0.6622333205008432
                b = 0.3377666794991567
            self.A_FeH = a * 10 ** aFe + b


    def __del__(self):
        """
        Clears a model instance.
        """

        del self.modes
        del self.glb
        del self.name

    def string_to_param(self, string, a=[]):
        """
        Return a parameter for an input string.

        :param string: string that indicates which parameter we're seeking
        :type string: string

        :param a: amplitude parameters which intervene in the surface correction
        :type a: array-like

        :return: the value of the parameter
        :rtype: float
        """

        if (string.startswith("log_")):
            return math.log10(self.string_to_param(string[4:]))
        if (string.startswith("ln_")):
            return math.log(self.string_to_param(string[3:]))
        if (string.startswith("exp_")):
            return math.exp(self.string_to_param(string[4:]))
        if (string == "Mass"):
            return self.glb[imass] / constants.solar_mass
        if (string == "Radius"):
            return self.glb[iradius] / constants.solar_radius
        if (string == "Luminosity"):
            return self.glb[iluminosity] / constants.solar_luminosity
        if (string == "Z"):
            return self.glb[iz0]
        if (string == "Y"):
            return 1.0 - self.glb[iz0] - self.glb[ix0]
        if (string == "X"):
            return self.glb[ix0]
        if (string == "Ys"):
            return 1.0 - self.glb[user_params_index["Zs"]] - self.glb[user_params_index["Xs"]]
        if (string == "Yc"):
            return 1.0 - self.glb[user_params_index["Zc"]] - self.glb[user_params_index["Xc"]]
        if (string == "zsx_s"):
            return self.zsx_s
        if (string == "zsx_0"):
            return self.zsx_0
        if (string == "Fe_H"):
            return self.FeH
        if (string == "Fe_H0"):
            return self.FeH0
        if (string == "M_H"):
            return self.MH
        if (string == "M_H0"):
            return self.MH0
        if (string == "dY_dZ"):
            return (1.0 - self.glb[iz0] - self.glb[ix0] - config.Yp) / self.glb[iz0]
        if (string == age_str):
            return self.glb[iage]
        if (string == age_adim_str):
            return self.glb[iage_adim]
        if (string == "Teff"):
            return self.glb[itemperature]
        if (string == "Dnu"):
            return self.find_large_separation() * self.glb[ifreq_ref]
        if (string == "Dnu_c"):
            return self.find_surface_corrected_large_separation(a) * self.glb[ifreq_ref]
        if (string == "numax"):
            return self.numax
        if (string == "Rho"):
            return 3.0 * self.glb[imass] / (4.0 * math.pi * self.glb[iradius] ** 3)
        if (string == "g"):
            return constants.G * self.glb[imass] / self.glb[iradius] ** 2
        if (string == "Theta_BCZ"):
            return self.glb[user_params_index["tau_BCZ"]] / (
                    2.0 * self.glb[user_params_index["tau"]])
        if (string == "Theta_He"):
            return self.glb[user_params_index["tau_He"]] / (
                    2.0 * self.glb[user_params_index["tau"]])
        if (string == "beta_Sonoi2015"):
            return self.beta_Sonoi2015
        if (string == "b_Kjeldsen2008"):
            return self.b_Kjeldsen2008

        # raises a KeyError if `string` isn't a valid key
        return self.glb[user_params_index[string]]

    def read_file(self, filename):
        """
        Read in a set of modes from a file.  This method will either call
        :py:meth:`read_file_CLES`, :py:meth:`read_file_MESA`, or
        :py:meth:`read_file_agsm` according to the value of the ``mode_format``
        variable in ``AIMS_configure.py``.

        :param filename: name of the file with the modes. The format of this file
                         is decided by the ``mode_format`` variable in
                         ``AIMS_configure.py``.

        :type filename: string

        :return: ``True`` if at least one frequency has been discarded (see note below).
        :rtype: boolean

        .. note::
          At this stage the frequencies should be expressed in :math:`\\mu\\mathrm{Hz}`.
          They will be non-dimensionalised in :py:func:`read_model_list_standard`.
        """

        if (config.mode_format == "simple"):  # for retrocompatibility
            return self.read_file_CLES(filename)
        if (config.mode_format == "CLES"):
            return self.read_file_CLES(filename)
        if (config.mode_format == "CLES_Mod"):
            return self.read_file_CLES_Mod(filename)
        if (config.mode_format == "MESA"):
            return self.read_file_MESA(filename)
        if (config.mode_format == "agsm"):
            return self.read_file_agsm(filename)
        if (config.mode_format == "PLATO"):
            return self.read_file_PLATO(filename)
        raise ValueError("Unrecognised format \"" + config.mode_format + "\".\n" +
                         "Please choose another value for mode_format in AIMS_configure.py")

    def read_file_CLES(self, filename):
        """
        Read in a set of modes from a file.  This uses the "simple" or "CLES"
        format as specified in the ``mode_format`` variable in ``AIMS_configure.py``.

        :param filename: name of the file with the modes.
          The file should contain a one-line header followed by five
          columns which correspond to l, n, frequency, unused, inertia.

        :type filename: string

        :return: ``True`` if at least one frequency has been discarded (see note below).
        :rtype: boolean

        .. note::
          - The fourth column is discarded.
          - Frequencies above `config.cutoff` * :math:`\\nu_{\\mathrm{cut-off}}` are discarded.
        """

        freqlim = config.cutoff * self.cutoff
        exceed_freqlim = False
        with open(filename) as freqfile:
            freqfile.readline()  # skip head
            mode_temp = []
            for line in freqfile:
                line = line.strip()
                columns = line.split()
                n = int(columns[1])
                freq = utilities.to_float(columns[2])
                # remove frequencies above AIMS_configure.cutoff*nu_{cut-off}
                if (freq > freqlim):
                    exceed_freqlim = True
                    continue
                if (config.npositive and (n < 0)):  # remove g-modes if need be
                    continue
                mode_temp.append((n, int(columns[0]), freq, utilities.to_float(columns[4])))

        self.modes = np.array(mode_temp, dtype=modetype)

        return exceed_freqlim

    def read_file_CLES_Mod(self, filename):
        """
        Read in a set of modes from a file.  This uses the "CLES_Mod" format as
        specified in the ``mode_format`` variable in ``AIMS_configure.py``.

        :param filename: name of the file with the modes.
          The file may contain a header (the header lines start with "#")
          followed by four columns which correspond to l, n, unused, frequency.

        :type filename: string

        :return: ``True`` if at least one frequency has been discarded (see note below).
        :rtype: boolean

        .. note::
          - The third column is discarded.
          - Frequencies above `config.cutoff` * :math:`\\nu_{\\mathrm{cut-off}}` are discarded.
        """

        freqlim = config.cutoff * self.cutoff
        exceed_freqlim = False
        with open(filename) as freqfile:
            mode_temp = []
            for line in freqfile:
                if line[0] == '#':  # Skip comments
                    continue
                line = line.strip()
                columns = line.split()
                n = int(columns[1])
                freq = utilities.to_float(columns[3]) * 10 ** 6  # freq tu muHz
                # remove frequencies above AIMS_configure.cutoff*nu_{cut-off}
                if (freq > freqlim):
                    exceed_freqlim = True
                    continue
                if (config.npositive and (n < 0)):  # remove g-modes if need be
                    continue
                mode_temp.append((n, int(columns[0]), freq, 1.0))  # No inertia in file

        self.modes = np.array(mode_temp, dtype=modetype)

        return exceed_freqlim

    def read_file_MESA(self, filename):
        """
        Read in a set of modes from a file.  This uses the GYRE (i.e. MESA) format as
        specified in the ``mode_format`` variable in ``AIMS_configure.py``.

        :param filename: name of the file with the modes.
          The file should contain a seven-line header followed by various
          columns which contain l, n, frequency, and inertia for each pulsation
          mode.

        :type filename: string

        :return: ``True`` if at least one frequency has been discarded (see note below).
        :rtype: boolean

        .. note::
          - Frequencies above `config.cutoff` * :math:`\\nu_{\\mathrm{cut-off}}` are discarded.
        """

        freqlim = config.cutoff * self.cutoff
        exceed_freqlim = False
        with open(filename) as freqfile:
            for i in range(6):
                freqfile.readline()
            mode_temp = []
            for line in freqfile:
                line = line.strip()
                columns = line.split()
                n = int(columns[1])  # n_pg
                # n = int(columns[2])  # n_p
                freq = utilities.to_float(columns[4])
                # remove frequencies above AIMS_configure.cutoff*nu_{cut-off}
                if (freq > freqlim):
                    exceed_freqlim = True
                    continue
                if (config.npositive and (n < 0)):  # remove g-modes if need be
                    continue
                mode_temp.append((n, int(columns[0]), freq, utilities.to_float(columns[7])))

        self.modes = np.array(mode_temp, dtype=modetype)

        return exceed_freqlim

    def read_file_agsm(self, filename):
        """
        Read in a set of modes from a file.  This uses the "agsm" format as
        specified in the ``mode_format`` variable in ``AIMS_configure.py``.

        :param filename: name of the file with the modes. This file is a
            binary fortran "agsm" file produced by the ADIPLS code.  See
            instructions to the ADIPLS code for a description of this
            format.

        :type filename: string

        :return: ``True`` if at least one frequency has been discarded (see note below).
        :rtype: boolean

        .. note::
          - Frequencies above `config.cutoff` * :math:`\\nu_{\\mathrm{cut-off}}` are discarded.
        """

        narr, larr, farr, iarr, nn, exceed_freqlim = \
            aims_fortran.read_file_agsm(filename, config.npositive, config.agsm_cutoff, \
                                        config.cutoff * self.cutoff)
        self.modes = np.array(list(zip(narr[0:nn], larr[0:nn], farr[0:nn], iarr[0:nn])), dtype=modetype)

        return exceed_freqlim

    def read_file_PLATO(self, filename):
        """
        Read in a set of modes from a file.  This uses the "PLATO" format as
        specified in the ``mode_format`` variable in ``AIMS_configure.py``.

        :param filename: name of the file with the modes.
          The file may contain a header (the header lines start with "#")
          followed by four columns which correspond to l, n, frequency, inertia.

        :type filename: string

        :return: ``False``
        :rtype: boolean

        .. note::
          - No mode selection is carried out, hence the return value
          - np.loadtxt is used to gain speed (the grid is quite large)
        """

        self.modes = np.loadtxt(filename, dtype=modetype, usecols=(1, 0, 2, 3))
        return False

    def read_modes_Aldo(self, norder, farray):
        """
        Read in a set of modes from a table.  This uses the Aldo format as
        specified in the ``mode_format`` variable in ``AIMS_configure.py``.

        :param norder: table which gives the starting n values and number of modes
                       per l value
        :type norder: 2D numpy int array

        :param farray: table with mode frequencies and inertias as follows:
          frequency1, inertia1, frequency2, inertia2 ...
        :type farray: 1D numpy float array

        :return: ``True`` if at least one frequency has been discarded (see note below).
        :rtype: boolean

        .. note::
          - Frequencies above `config.cutoff` * :math:`\\nu_{\\mathrm{cut-off}}` are discarded.
        """

        freqlim = config.cutoff * self.cutoff
        exceed_freqlim = False
        mode_temp = []
        ii = 0
        for l in range(4):
            for i in range(norder[l, 1]):
                n = i + norder[l, 0]
                freq = farray[ii]

                # remove zero freuqency modes (these modes are not calculated):
                if (abs(freq) < 1e-10):
                    ii += 2
                    continue

                # remove frequencies above AIMS_configure.cutoff*nu_{cut-off}
                if (freq > freqlim):
                    exceed_freqlim = True
                    ii += 2
                    continue

                # remove gravity modes if need be
                if (config.npositive and (n < 0)):
                    ii += 2
                    continue

                # we can add this mode:
                mode_temp.append((n, l, freq, farray[ii + 1]))
                self.modes = np.array(mode_temp, dtype=modetype)
                ii += 2

        return exceed_freqlim

    def write_file_simple(self, filename):
        """
        Write a set of modes into a file using the "simple" format as
        described in :py:meth:`read_file_simple`.

        :param filename: name of the file where the modes should
                         be written.

        :type filename: string

        .. note::
          - The output frequencies are expressed in :math:`\\mathrm{\\mu Hz}`
        """

        with open(filename, "w") as output:
            # write header
            output.write("# %1s %3s %22s %6s %22s\n" % ("l", "n", "nu_theo (muHz)", "unused", "Inertia"))
            for i in range(self.modes.shape[0]):
                output.write("  %1d %3d %22.15e    0.0 %22.15e\n" % (
                    self.modes["l"][i],
                    self.modes["n"][i],
                    self.modes["freq"][i] * self.glb[ifreq_ref],
                    self.modes["inertia"][i]))

    def append_modes(self, modes):
        """
        Append a list of modes to the model.

        :param modes: list of modes which are in the form of tuples: (n,l,freq,inertia).
        :type modes: list of (int, int, float, float)
        """
        self.modes = np.append(self.modes, np.array(modes, dtype=modetype))

    def sort_modes(self):
        """
        Sort the modes by l, then n, then freq.
        """
        # sorts by l, then n, then freq
        ind = np.lexsort((self.modes['freq'], self.modes['n'], self.modes['l']))
        self.modes = np.array([self.modes[i] for i in ind], dtype=modetype)

    def remove_duplicate_modes(self):
        """
        Remove duplicate modes.

        Modes are considered to be duplicate if they have the same l
        and n values (regardless of frequency).

        :return: ``True`` if at least one mode has been removed.
        :rtype: boolean

        .. warning::
          This method assumes the modes are sorted.
        """

        ind = []
        for i in range(len(self.modes) - 1, 1, -1):
            if ((self.modes['n'][i] == self.modes['n'][i - 1]) and \
                    (self.modes['l'][i] == self.modes['l'][i - 1])):
                ind.append(i)
        self.modes = np.delete(self.modes, ind)
        return (len(ind) == 0)

    def get_age(self):
        """
        Return age of stellar model.

        This is useful for sorting purposes.

        :return: the age of the model
        :rtype: float
        """
        return self.glb[iage]

    def get_freq(self, surface_option=None, a=[]):
        """
        Obtain model frequencies, with optional frequency corrections.

        :param surface_option: specifies the type of surface correction.  Options include:

          - ``None``: no corrections are applied
          - ``"Kjeldsen2008"``: apply a correction based on Kjeldsen et al. (2008)
          - ``"Kjeldsen2008_scaling"``: apply a correction based on Kjeldsen et al. (2008).
                                        The exponent is based on a scaling relation from
                                        Sonoi et al. (2015).
          - ``"Kjeldsen2008_2"``: apply a correction based on Kjeldsen et al. (2008).
                                        The exponent is a free parameter.
          - ``"Ball2014"``:     apply a one-term correction based on Ball and Gizon (2014)
          - ``"Ball2014_2"``:   apply a two-term correction based on Ball and Gizon (2014)
          - ``"Sonoi2015"``: apply a correction based on Sonoi et al. (2015)
          - ``"Sonoi2015_scaling"``: apply a correction based on Sonoi et al. (2015)
                                        The exponent is based on a scaling relation from
                                        Sonoi et al. (2015).
          - ``"Sonoi2015_2"``: apply a correction based on Sonoi et al. (2015)
                                        The exponent is a free parameter.

        :param a: amplitude parameters which intervene in the surface correction

        :type surface_option: string
        :type a: array-like

        :return: models frequencies (including surface corrections)
        :rtype: np.array

        .. note::
          If surface_option==None or a==[], the original frequencies are
          returned (hence modifying them modifies the :py:class:`Model` object).
        """

        if (surface_option is None) or (len(a) == 0):
            return self.modes['freq']
        return self.modes['freq'] + self.get_surface_correction(surface_option, a)

    def get_surface_correction(self, surface_option, a):
        """
        Obtain corrections on model frequencies (these corrections should be
        *added* to the *theorectical* frequencies).

        :param surface_option: specifies the type of surface correction.
          Options include:

          - ``None``: no corrections are applied
          - ``"Kjeldsen2008"``: apply a correction based on Kjeldsen et al. (2008)
          - ``"Kjeldsen2008_scaling"``: apply a correction based on Kjeldsen et al. (2008).
                                        The exponent is based on a scaling relation from
                                        Sonoi et al. (2015).
          - ``"Kjeldsen2008_2"``: apply a correction based on Kjeldsen et al. (2008).
                                        The exponent is a free parameter.
          - ``"Ball2014"``:     apply a one-term correction based on Ball and Gizon (2014)
          - ``"Ball2014_2"``:   apply a two-term correction based on Ball and Gizon (2014)
          - ``"Sonoi2015"``: apply a correction based on Sonoi et al. (2015)
          - ``"Sonoi2015_scaling"``: apply a correction based on Sonoi et al. (2015)
                                        The exponent is based on a scaling relation from
                                        Sonoi et al. (2015).
          - ``"Sonoi2015_2"``: apply a correction based on Sonoi et al. (2015)
                                        The exponent is a free parameter.

        :param a: parameters which intervene in the surface correction.  According to
          the correction they take on the following meanings:

          - ``"Kjeldsen2008"``:         a[0]*freq**b_Kjeldsen2008
          - ``"Kjeldsen2008_scaling"``: a[0]*freq**b_scaling
          - ``"Kjeldsen2008_2"``:       a[0]*freq**a[1]
          - ``"Ball2014"``:             a[0]*freq**3/I
          - ``"Ball2014_2"``:           a[0]*freq**3/I + a[1]/(freq*I)
          - ``"Sonoi2015"``:            a[0]*[1 - 1/(1 + (nu/numax)**beta_Sonoi2015)]
          - ``"Sonoi2015_scaling"``:    a[0]*[1 - 1/(1 + (nu/numax)**beta_scaling)]
          - ``"Sonoi2015_2"``:          a[0]*[1 - 1/(1 + (nu/numax)**a[1])]

        :type surface_option: string
        :type a: array-like

        :return: surface corrections on the model frequencies
        :rtype: np.array

        .. note::
          The array operations lead to the creation of a new array with the
          result, which avoids modifications of the original frequencies and inertias.
        """

        # easy exit:
        if (surface_option is None):
            return np.zeros(len(self.modes), dtype=ftype)

        freq = self.glb[ifreq_ref] * self.modes['freq'] / self.numax
        if (surface_option == "Kjeldsen2008"):
            return a[0] * freq ** config.b_Kjeldsen2008
        if (surface_option == "Kjeldsen2008_scaling"):
            return a[0] * freq ** self.b_Kjeldsen2008
        if (surface_option == "Kjeldsen2008_2"):
            return a[0] * freq ** a[1]
        if (surface_option == "Ball2014"):
            return a[0] * freq ** 3 / self.modes['inertia']
        if (surface_option == "Ball2014_2"):
            return a[0] * freq ** 3 / self.modes['inertia'] \
                + a[1] / (freq * self.modes['inertia'])
        if (surface_option == "Sonoi2015"):
            return a[0] * (1.0 - 1.0 / (1.0 + freq ** config.beta_Sonoi2015))
        if (surface_option == "Sonoi2015_scaling"):
            return a[0] * (1.0 - 1.0 / (1.0 + freq ** self.beta_Sonoi2015))
        if (surface_option == "Sonoi2015_2"):
            return a[0] * (1.0 - 1.0 / (1.0 + freq ** a[1]))
        raise ValueError("Unknown surface correction: " + surface_option)

    @property
    def b_Kjeldsen2008(self):
        """
        Return the exponent for the Kjeldsen et al. (2008) surface correction
        recipe, as calculated based on the Sonoi et al. (2015) scaling relation.

        :return: the Kjeldsen et al. exponent
        :rtype: float
        """
        return 10.0 ** (-3.16 * self.string_to_param("log_Teff") + 0.184 * self.string_to_param("log_g") + 11.7)

    @property
    def beta_Sonoi2015(self):
        """
        Return the exponent for the Sonoi et al. (2015) surface correction
        recipe, as calculated based on the Sonoi et al. (2015) scaling relation.

        :return: the Kjeldsen et al. exponent
        :rtype: float
        """
        return 10.0 ** (-3.86 * self.string_to_param("log_Teff") + 0.235 * self.string_to_param("log_g") + 14.2)

    def multiply_modes(self, constant):
        """
        Multiply the frequencies by constant.

        :param constant: constant by which the mode frequencies are multiplied
        :type constant: float
        """

        # NOTE: inertias are non-dimensionless, so they don't change.
        # (this has been tested numerically)
        self.modes['freq'] *= constant

    def find_mode(self, ntarget, ltarget):
        """
        Find a mode with specific n and l values.

        :param ntarget: target n value
        :param ltarget: target l value
        :type ntarget: int
        :type ltarget: int

        :return: the frequency of the mode
        :rtype: float
        """

        size = len(self.modes)
        for i in range(size):
            if ((self.modes['n'][i] == ntarget) and (self.modes['l'][i] == ltarget)):
                return self.modes['freq'][i]
        return np.nan

    def find_mode_range(self):
        """
        Find n and l ranges of the modes in the model.

        :return: the n and l ranges of the modes
        :rtype: int, int, int, int
        """

        if (len(self.modes) < 1):
            return -1, -1, -1, -1
        nmin = np.nanmin(self.modes['n'])
        nmax = np.nanmax(self.modes['n'])
        lmin = np.nanmin(self.modes['l'])
        lmax = np.nanmax(self.modes['l'])
        return nmin, nmax, lmin, lmax

    def find_large_separation(self):
        """
        Find large frequency separation using only radial modes.

        :return: the large frequency separation
        :rtype: float
        """

        ind = (self.modes['l'] == 0).flat
        n = self.modes['n'].compress(ind)
        if (len(np.unique(n)) < 2):
            return np.nan
        freq = self.modes['freq'].compress(ind)
        # Use numax in uHz to calculate sigma, then convert to dimless
        numax_dimless = self.numax / self.glb[ifreq_ref]
        sigma_dimless = 0.66 * self.numax ** 0.88 / self.glb[ifreq_ref]
        weights = np.exp(-((freq - numax_dimless) / sigma_dimless) ** 2)
        coeff = np.polyfit(n, freq, deg=1, w=weights)
        return coeff[0]

    def find_surface_corrected_large_separation(self, a=[]):
        """
        Find large frequency separation using only radial modes and surface corrections.

        :param a: Surface correction parameters.

        :return: the large frequency separation
        :rtype: float
        """

        ind = (self.modes['l'] == 0).flat
        n = self.modes['n'].compress(ind)
        if (len(np.unique(n)) < 2):
            return np.nan

        # Weight large sep by envelope
        freq = self.get_freq(config.surface_option, a).compress(ind)
        numax_dimless = self.numax / self.glb[ifreq_ref]
        sigma_dimless = 0.66 * self.numax ** 0.88 / self.glb[ifreq_ref]
        weights = np.exp(-((freq - numax_dimless) / sigma_dimless)**2)
        coeff = np.polyfit(n, freq, deg=1, w=weights)
        return coeff[0]

    def find_epsilon(self, ltarget):
        """
        Find epsilon, the constant offset in a simplified version of Tassoul's
        asymptotic formula:

        :math:`\\nu_n = \\Delta \\nu (n + \\varepsilon)`

        :param ltarget: target l value.  Only modes with this l value will be
          used in obtaining epsilon.
        :type ltarget: int

        :return: the constant offset
        :rtype: float
        """

        dnu = self.find_large_separation()
        one = n = nu = 0.0
        for i in range(len(self.modes)):
            if (self.modes['l'][i] != ltarget):
                continue
            one += 1.0
            n += self.modes['n'][i]
            nu += self.modes['freq'][i]
        if (one == 0.0):
            return 0.0
        else:
            return (nu / dnu - n) / one

    @property
    def FeH(self):
        """
        Find [Fe/H] value for model.

        The conversion from (Xs,Zs) to [Fe/H] is performed using the
        following formula:

            :math:`\\mathrm{[Fe/H] = \\frac{[M/H]}{A_{FeH}}  \
                                   = \\frac{1}{A_{FeH}} \\log_{10} \
                                     \\left(\\frac{z/x}{z_{\\odot}/x_{\\odot}} \\right)}`

        :return: the :math:`\\mathrm{[Fe/H]}` value
        :rtype: float

        .. note::
          The relevant values are given in :py:mod:`constants`
        """
        try:
            return self.MH / self.A_FeH
        except ValueError:
            return log0  # a rather low value

    @property
    def FeH0(self):
        """
        Find initial [Fe/H] value for model.

        The conversion from (X,Z) to [Fe/H] is performed using the
        following formula:

            :math:`\\mathrm{[Fe/H] = \\frac{[M/H]}{A_{FeH}}  \
                                   = \\frac{1}{A_{FeH}} \\log_{10} \
                                     \\left(\\frac{z/x}{z_{\\odot}/x_{\\odot}} \\right)}`

        :return: the initial :math:`\\mathrm{[Fe/H]}` value
        :rtype: float

        .. note::
          The relevant values are given in :py:mod:`constants`
        """

        try:
            return self.MH0 / self.A_FeH
        except ValueError:
            return log0  # a rather low value

    @property
    def MH(self):
        """
        Find [M/H] value for model.

        The conversion from (Xs,Zs) to [M/H] is performed using the
        following formula:

            :math:`\\mathrm{[M/H] = \\log_{10} \\left(\\frac{z/x}{z_{\\odot}/x_{\\odot}} \\right)}`

        :return: the :math:`\\mathrm{[M/H]}` value
        :rtype: float

        .. note::
          The relevant values are given in :py:mod:`constants`
        """
        available_keys = user_params_index.keys()
        have_Zs = False
        have_Ys = False
        have_Xs = False
        have_ZXs = False
        if 'Zs' in available_keys:
            Zs = self.glb[user_params_index["Zs"]]
            have_Zs = True
        if 'Ys' in available_keys:
            Ys = self.glb[user_params_index["Ys"]]
            have_Ys = True
        if 'Xs' in available_keys:
            Xs = self.glb[user_params_index["Xs"]]
            have_Xs = True
        if 'ZXs' in available_keys:
            ZXs = self.glb[user_params_index["ZXs"]]
            have_ZXs = True
        if 'ZHsurf' in available_keys:
            ZXs = self.glb[user_params_index["ZHsurf"]]
            have_ZXs = True

        if not have_Zs:
            if have_Xs and have_Ys:
                Zs = 1 - self.glb[user_params_index["Xs"]] - self.glb[user_params_index["Ys"]]
                have_Zs = True
            if have_ZXs and have_Xs:
                Zs = ZXs * Xs
                have_Zs = True

        try:
            return math.log10(Zs * config.solar_x / (Xs * config.solar_z))
        except ValueError:
            return log0  # a rather low value

    @property
    def MH0(self):
        """
        Find initial [M/H] value for model.

        The conversion from (X,Z) to [M/H] is performed using the
        following formula:

            :math:`\\mathrm{[M/H] = \\log_{10} \\left(\\frac{z/x}{z_{\\odot}/x_{\\odot}} \\right)}`

        :return: the initial :math:`\\mathrm{[M/H]}` value
        :rtype: float

        .. note::
          The relevant values are given in :py:mod:`constants`
        """

        try:
            return math.log10(self.glb[iz0] * config.solar_x / (self.glb[ix0] * config.solar_z))
        except ValueError:
            return log0  # a rather low value

    @property
    def zsx_s(self):
        """
        Find the Zs/Xs value

        :return: the Zs/Xs value
        :rtype: float
        """

        return self.glb[user_params_index["Zs"]] / self.glb[user_params_index["Xs"]]

    @property
    def zsx_0(self):
        """
        Find the Z0/X0 value

        :return: the Z0/X0 value
        :rtype: float
        """
        return self.glb[iz0] / self.glb[ix0]

    @property
    def numax(self):
        """
        Find :math:`\\nu_{\\mathrm{max}}` for model.

        The :math:`\\nu_{\\mathrm{max}}` value is obtained from the
        following scaling relation:

            :math:`\\frac{\\nu_{\\mathrm{max}}}{\\nu_{\\mathrm{max},\\odot}} \
                      = \\left(\\frac{M}{M_{\\odot}}\\right)                \
                        \\left(\\frac{R}{R_{\\odot}}\\right)^{-2}           \
                        \\left(\\frac{T_{\\mathrm{eff}}}{T_{\\mathrm{eff},\\odot}}\\right)^{-1/2}`

        :return: the :math:`\\nu_{\\mathrm{max}}` value
        :rtype: float

        .. note::
          The relevant values are given in :py:mod:`constants`
        """

        return constants.solar_numax * (self.glb[imass] / constants.solar_mass) \
            / ((self.glb[iradius] / constants.solar_radius) ** 2 \
               * math.sqrt(self.glb[itemperature] / constants.solar_temperature))

    @property
    def cutoff(self):
        """
        Find :math:`\\nu_{\\mathrm{cut-off}}` for model.

        The :math:`\\nu_{\\mathrm{cut-off}}` value is obtained from the
        following scaling relation:

            :math:`\\frac{\\nu_{\\mathrm{cut-off}}}{\\nu_{\\mathrm{cut-off},\\odot}} \
                      = \\left(\\frac{M}{M_{\\odot}}\\right)                        \
                        \\left(\\frac{R}{R_{\\odot}}\\right)^2                      \
                        \\left(\\frac{T_{\\mathrm{eff}}}{T_{\\mathrm{eff},\\odot}}\\right)^{-1/2}`

        :return: the :math:`\\nu_{\\mathrm{cut-off}}` value
        :rtype: float

        .. note::
          The relevant values are given in :py:mod:`constants`
        """

        return constants.solar_cutoff * (self.glb[imass] / constants.solar_mass) \
            / ((self.glb[iradius] / constants.solar_radius) ** 2 \
               * math.sqrt(self.glb[itemperature] / constants.solar_temperature))

    def freq_sorted(self):
        """
        Check to see if the frequencies are in ascending order for each l value.

        :return: ``True`` if the frequencies are in ascending order.
        :rtype: boolean
        """

        for i in range(len(self.modes) - 1):
            # if (self.modes['l'][i] > 0): continue
            if (self.modes['l'][i] != self.modes['l'][i + 1]):
                continue
            if (self.modes['freq'][i] > self.modes['freq'][i + 1]):
                return False
        return True

    def print_me(self):
        """ Print classical and seismic characteristics of the model to standard output."""

        print("----- Model: " + str(self.name) + " -----")
        print("Mass (in M_sun):              %.5f" % (self.glb[imass] / constants.solar_mass))
        print("Radius (in R_sun):            %.5f" % (self.glb[iradius] / constants.solar_radius))
        print("Reference frequency (in uHz): %.3f" % self.glb[ifreq_ref])
        print("Temperature (in K):           %.1f" % self.glb[itemperature])
        print("Luminosity (in L_sun):        %.3g" % (self.glb[iluminosity] / constants.solar_luminosity))
        print("Age (in Myrs):                %.2f" % self.glb[iage])
        print("Age parameter:                %.2f" % self.glb[iage_adim])
        print("Z:                            %.4f" % self.glb[iz0])
        print("X:                            %.4f" % self.glb[ix0])
        for (name, latex_name) in config.user_params:
            print("{0:29} {1:.5e}".format(name, self.glb[user_params_index[name]]))
        print("Modes (in muHz):")
        size = self.modes.shape[0]
        for i in range(size):
            print("  (n,l,freq,IK) = (%d, %d, %.15f, %.5e)" % \
                  (self.modes['n'][i], self.modes['l'][i], \
                   self.modes['freq'][i] * self.glb[ifreq_ref], \
                   self.modes['inertia'][i]))


class Track:
    """
    An evolutionary track.
    """

    def __init__(self, aModel, grid_params):
        """
        :param aModel: first model to be added to evolutionary track (it does
          not need to be the youngest model in an evolutionary sequence).  This
          Model is used to obtain the relevant parameters for the evolutionary
          track (as given by the :py:data:`grid_params` variable).

        :param grid_params: list of strings which are the names of the
          parameters which describe the evolutionary track.

        :type aModel: :py:class:`Model`
        :type grid_params: list of strings
        """

        self.grid_params = grid_params
        """Names of the parameters used to construct the grid"""

        self.params = np.array(utilities.my_map(aModel.string_to_param, self.grid_params))
        """Parameters which characterise this evolutionary track"""

        self.nmodes = len(aModel.modes)
        """Total number pulsation modes from all of the models in this evolutionary track"""

        self.names = [aModel.name, ]
        """ List of model names """

        self.glb = aModel.glb.reshape((1, nglb))
        """ Global properties of the models """

        self.modes = aModel.modes
        """ array containing the modes (n, l, freq, inertia) of all of the models """

        self.mode_indices = [0, len(aModel.modes)]
        """ starting indices in :py:data:`Track.modes` array corresponding to a given model """

        del aModel  # clean up to avoid using too much memory

    def __del__(self):
        """
        Remove a track and clear the arrays it contained.
        """

        del self.grid_params
        del self.params
        del self.nmodes
        del self.names
        del self.glb
        del self.modes
        del self.mode_indices

    def append(self, aModel):
        """
        Append a model to the evolutionary track.

        :param aModel: model which is being appended to the track
        :type aModel: :py:class:`Model`
        """

        self.nmodes += len(aModel.modes)
        self.names.append(aModel.name)
        self.glb = np.append(self.glb, aModel.glb.reshape((1, nglb)), axis=0)
        self.modes = np.append(self.modes, aModel.modes)
        self.mode_indices.append(self.nmodes)
        del aModel  # clean up to avoid using too much memory

    def append_track(self, aTrack):
        """
        Append a track to the current evolutionary track (i.e. combine the two tracks),
        and remove the track which has been appended.  The resultant track is then
        sorted according to age.

        :param aModel: model which is being appended to the track
        :type aModel: :py:class:`Model`
        """

        self.names += aTrack.names
        self.glb = np.concatenate((self.glb, aTrack.glb), axis=0)
        self.modes = np.concatenate((self.modes, aTrack.modes), axis=0)
        self.mode_indices += [self.nmodes + i for i in aTrack.mode_indices[1:]]
        self.nmodes += aTrack.nmodes
        self.sort()

    def matches(self, aModel, params_tol=None):
        """
        Check to see if a model matches the evolutionary track and can therefore
        be included in the track.

        :param aModel: input model being tested
        :type aModel: :py:class:`Model`

        :param params_tol: the tolerance on each parameter. ``None`` means no
                 tolerance (i.e. the parameters have to be exactly the same).
        :type params_tol: float np.array

        :return: ``True`` only if all of the parameters of the input model match
                 those of the evolutionary track within the provided tolerances
        :rtype: boolean
        """

        params_bis = np.array(utilities.my_map(aModel.string_to_param, self.grid_params))
        if (params_tol is None):
            return np.all(self.params == params_bis)
        else:
            return np.all(np.abs(self.params - params_bis) <= params_tol)

    def sort(self):
        """Sort models within evolutionary track according to age."""

        index = np.argsort(self.glb[:, iage])

        # sort model names
        temp = utilities.my_map(lambda i: self.names[i], index)
        del self.names
        self.names = temp

        # sort global parameters
        temp = self.glb[index]
        del self.glb
        self.glb = temp

        # sort modes and mode indices:
        temp = np.empty(self.modes.shape, dtype=modetype)
        ntot = 0
        temp2 = [ntot, ]
        for i in index:
            nmodes = self.mode_indices[i + 1] - self.mode_indices[i]
            temp[ntot:ntot + nmodes] = self.modes[self.mode_indices[i]:self.mode_indices[i + 1]]
            ntot += nmodes
            temp2.append(ntot)
        del self.modes
        del self.mode_indices
        self.modes = temp
        self.mode_indices = temp2

    def is_sorted(self):
        """
        Check to see of models are in ascending order according to age.

        :return: ``True`` if the models ar in order of increasing age
        :rtype: boolean
        """

        return all(self.glb[i, iage] < self.glb[i + 1, iage] for i in range(len(self.names) - 1))

    def is_sorted_adim(self):
        """
        Check to see of models are in ascending order according to age.

        :return: ``True`` if the models ar in order of increasing age
        :rtype: boolean
        """

        return all(self.glb[i, iage_adim] < self.glb[i + 1, iage_adim] for i in range(len(self.names) - 1))

    def freq_sorted(self, imodel):
        """
        Check to see if the frequencies of model i are in ascending order
        for each l value.

        :param imodel: the number of the model to be checked
        :type imodel: int

        :return: ``True`` if the frequencies are in ascending order.
        :rtype: boolean
        """

        for i in range(self.mode_indices[imodel], self.mode_indices[imodel + 1] - 1):
            # if (self.modes['l'][i] > 0): continue
            if (self.modes['l'][i] != self.modes['l'][i + 1]):
                continue
            if (self.modes['freq'][i] > self.modes['freq'][i + 1]):
                return False
        return True

    def find_epsilon(self, imodel, ltarget):
        """
        Find epsilon, the constant offset in a simplified version of Tassoul's
        asymptotic formula:

        :math:`\\nu_n = \\Delta \\nu (n + \\varepsilon)`

        :param imodel: model number
        :type imodel: ltarget

        :param ltarget: target l value.  Only modes with this l value will be
          used in obtaining epsilon.
        :type ltarget: int

        :return: the constant offset
        :rtype: float
        """

        dnu = self.find_large_separation(imodel)
        one = n = nu = 0.0
        for i in range(self.mode_indices[imodel], self.mode_indices[imodel + 1]):
            if (self.modes['l'][i] != ltarget):
                continue
            one += 1.0
            n += self.modes['n'][i]
            nu += self.modes['freq'][i]
        if (one == 0.0):
            return 0.0
        else:
            return (nu / dnu - n) / one

    def find_large_separation(self, imodel):
        """
        Find large frequency separation using only radial modes.

        :param imodel: model number
        :type imodel: ltarget

        :return: the large frequency separation
        :rtype: float
        """

        istart = self.mode_indices[imodel]
        istop = self.mode_indices[imodel + 1]
        ind = (self.modes['l'][istart:istop] == 0).flat
        n = self.modes['n'][istart:istop].compress(ind)
        if (len(n) < 2):
            return np.nan
        freq = self.modes['freq'][istart:istop].compress(ind)
        numax = (constants.solar_numax * (self.glb[imodel, imass] / constants.solar_mass) /
                 ((self.glb[imodel, iradius] / constants.solar_radius) ** 2 *
                  np.sqrt(self.glb[imodel, itemperature] / constants.solar_temperature)))
        numax_dimless = numax / self.glb[imodel, ifreq_ref]
        sigma_dimless = 0.66 * numax ** 0.88 / self.glb[imodel, ifreq_ref] # Use numax in uHz, then convert to dimless
        weights = np.exp(-((freq - numax_dimless) / sigma_dimless) ** 2)
        coeff = np.polyfit(n, freq, deg=1, w=weights)

        return coeff[0]

    def duplicate_ages(self):
        """
        Check to see if the evolutionary track contains models with duplicate ages.

        :return: ``True`` if there are duplicate age(s)
        :rtype: boolean

        .. warning::
            This method should only be applied after the track has been
            sorted.
        """

        return any(self.glb[i, iage] == self.glb[i + 1, iage] for i in range(len(self.names) - 1))

    def remove_duplicate_ages(self):
        """
        Removes models with duplicate ages from the evolutionary track.

        :return: ``True`` if models with duplicate age(s) have been removed
        :rtype: boolean

        .. warning::
            This method should only be applied after the track has been
            sorted.
        """

        # produce list of models to be removed:
        remove_models = [i for i in range(1, len(self.names)) \
                         if self.glb[i - 1, iage] == self.glb[i, iage]]

        return self.remove_models(remove_models)


    def remove_models(self, remove_models):
        if (len(remove_models) == 0):
            return False
        else:
            # remove duplicate models:
            mode_indices_rm = []
            for i in remove_models:
                mode_indices_rm += list(range(self.mode_indices[i], self.mode_indices[i + 1]))

            for i in remove_models[::-1]:
                del self.names[i]
            self.glb = np.delete(self.glb, remove_models, 0)
            self.modes = np.delete(self.modes, mode_indices_rm, 0)

            ntot = self.mode_indices[1]
            new_indices = [0, ntot]
            # NOTE: the first model is kept (by construction)
            for i in range(1, len(self.mode_indices) - 1):
                if (i not in remove_models):
                    ntot += self.mode_indices[i + 1] - self.mode_indices[i]
                    new_indices.append(ntot)
            self.mode_indices = new_indices
            self.nmodes = ntot

            return True

    def age_adim_to_age(self):
        """
        Replace dimensionless age parameter by the physical age.

        .. warning::
            This method should only be applied after the track has been
            sorted (according to age).
        """
        self.glb[:, iage_adim] = self.glb[:, iage]

    def age_adim_to_scale_age(self):
        """
        Replace dimensionless age parameter by the physical age scaled from 0 to 1.

        .. warning::
            This method should only be applied after the track has been
            sorted (according to age).
        """
        age_start = self.glb[0, iage]
        age_stop = self.glb[-1, iage]
        self.glb[:, iage_adim] = (self.glb[:, iage] - age_start) / (age_stop - age_start)

    def age_adim_to_scale_Xc(self):
        """
        Replace dimensionless age parameter by the central hydrogen abundance, Xc,
        scaled from 0 to 1.

        .. warning::
            This method should only be applied after the track has been
            sorted (according to age).
        """

        # Raises a key error if Xc isn't in the grid,
        # which should be informative enough
        ixc = user_params_index["Xc"]

        Xc_start = self.glb[0, ixc]
        Xc_stop = self.glb[-1, ixc]
        self.glb[:, iage_adim] = (self.glb[:, ixc] - Xc_start) / (Xc_stop - Xc_start)

    def interpolate_model(self, age):
        """
        Return a model at a given age which is obtained using linear interpolation.

        :param age: age of desired model in :math:`\\mathrm{Myrs}`
        :type age: float

        :return: the interpolated model
        :rtype: :py:class:`Model`

        .. warning::
          This method assumes the track is sorted, since it applies
          a binary search algorithm for increased efficiency.
        """

        # easy exit:
        if (age < self.glb[0, iage]):
            return None
        if (age > self.glb[-1, iage]):
            return None

        istart = 0
        istop = len(self.names) - 1
        while (istop > istart + 1):
            itemp = (istop + istart) // 2
            if (age < self.glb[itemp, iage]):
                istop = itemp
            else:
                istart = itemp
        mu = (age - self.glb[istart, iage]) \
             / (self.glb[istop, iage] - self.glb[istart, iage])
        aFe = get_aFe(self)
        return combine_models(
            Model(self.glb[istart], _modes=self.modes[self.mode_indices[istart]:self.mode_indices[istart + 1]], aFe=aFe),
            1.0 - mu, \
            Model(self.glb[istop], _modes=self.modes[self.mode_indices[istop]:self.mode_indices[istop + 1]], aFe=aFe), mu)

    def find_combination(self, age, coef):
        """
        Return a model combination at a given age which is obtained using linear interpolation.

        :param age: age of desired model in :math:`\\mathrm{Myrs}`
        :param coef: coefficient which multiplies this combination

        :type age: float
        :type coef: float

        :return: pairs composed of an interpolation coefficient and a model name
        :rtype: tuple of (float, string)

        .. warning::
          This method assumes the track is sorted, since it applies
          a binary search algorithm for increased efficiency.
        """

        # easy exit:
        if (age < self.glb[0, iage]):
            return None
        if (age > self.glb[-1, iage]):
            return None

        istart = 0
        istop = len(self.names) - 1
        while (istop > istart + 1):
            itemp = (istop + istart) // 2
            if (age < self.glb[itemp, iage]):
                istop = itemp
            else:
                istart = itemp
        mu = (age - self.glb[istart, iage]) \
             / (self.glb[istop, iage] - self.glb[istart, iage])
        return ((coef * (1.0 - mu), self.names[istart]), (coef * mu, self.names[istop]))

    def find_modes(self, ntarget, ltarget):
        """
        Return two lists, one with the ages of the models and the other
        with the mode non-dimensional frequencies corresponding to target n and l values.

        This function is useful for seeing how the frequency of a particular
        mode changes with stellar age.

        :param ntarget: target n value
        :param ltarget: target l value

        :type ntarget: int
        :type ltarget: int

        :return: lists of ages and frequencies
        :rtype: list, list
        """

        ages = list(self.glb[:, iage])
        freqs = []
        for i in range(len(self.names)):
            for j in range(self.mode_indices[i], self.mode_indices[i + 1]):
                if ((self.modes['n'][j] == ntarget) and (self.modes['l'][j] == ltarget)):
                    freqs.append(self.modes['freq'][j])
            else:
                freqs.append(np.nan)

        return ages, freqs

    def find_modes_dim(self, ntarget, ltarget):
        """
        Return two lists, one with the ages of the models and the other
        with the mode dimensional frequencies corresponding to target n and l values.

        This function is useful for seeing how the frequency of a particular
        mode changes with stellar age.

        :param ntarget: target n value
        :param ltarget: target l value

        :type ntarget: int
        :type ltarget: int

        :return: lists of ages and frequencies
        :rtype: list, list
        """

        ages = list(self.glb[:, iage])
        freqs = []
        for i in range(len(self.names)):
            for j in range(self.mode_indices[i], self.mode_indices[i + 1]):
                if ((self.modes['n'][j] == ntarget) and (self.modes['l'][j] == ltarget)):
                    freqs.append(self.modes['freq'][j] * self.glb[i, ifreq_ref])
            else:
                freqs.append(np.nan)

        return ages, freqs

    def find_mode_range(self):
        """
        Find n and l ranges of modes in models

        :return: the n and l ranges
        :rtype: int, int, int, int
        """

        if (len(self.modes) < 1):
            return -1, -1, -1, -1
        else:
            return np.min(self.modes['n']), np.max(self.modes['n']), \
                np.min(self.modes['l']), np.max(self.modes['l'])

    @property
    def age_range(self):
        """
        Provides the age range for an evolutionary track.
        """

        return abs(self.glb[-1, iage] - self.glb[0, iage])

    @property
    def age_lower(self):
        """
        Provides the lowest age an evolutionary track.
        """

        return min(self.glb[0, iage], self.glb[-1, iage])

    @property
    def age_upper(self):
        """
        Provides the highest age for an evolutionary track.
        """

        return max(self.glb[0, iage], self.glb[-1, iage])

    def test_interpolation(self, nincr):
        """
        Test accuracy of interpolation along evolutionary track.

        This method removes every other model and retrieves its frequencies
        by interpolation from neighbouring models.  The accuracy of the
        interpolated frequencies and global parameters are tested by carrying
        out comparisons with the original models.

        :param nincr: increment with which to carry out the interpolation.
          By comparing results for different values of ``nincr``, one can
          evaluate how the interpolation error depends on the size of the
          interval over which the interpolation is carried out.
        :type nincr: int

        :return: the interpolation errors
        :rtype: np.array
        """

        # initialisation
        nmodels = len(self.names)
        ndim = len(self.params) + 1
        result = np.zeros((nmodels - 2 * nincr, ndim + nglb + 6), dtype=gtype)

        aFe = get_aFe(self)
        # loop through all models:
        for i in range(nincr, nmodels - nincr):
            # carry out interpolation
            mu = (self.glb[i, iage] - self.glb[i - nincr, iage]) \
                 / (self.glb[i + nincr, iage] - self.glb[i - nincr, iage])
            aModel = combine_models(Model(self.glb[i - nincr], aFe=aFe, _modes=self.modes[
                                                                      self.mode_indices[i - nincr]:self.mode_indices[
                                                                          i - nincr + 1]]), 1.0 - mu, \
                                    Model(self.glb[i + nincr], aFe=aFe, _modes=self.modes[
                                                                      self.mode_indices[i + nincr]:self.mode_indices[
                                                                          i + nincr + 1]]), mu)

            result[i - nincr, 0:ndim - 1] = self.params
            result[i - nincr, ndim - 1] = self.glb[i, iage]
            result[i - nincr, ndim:ndim + nglb + 6] = compare_models(aModel, Model(self.glb[i], aFe=aFe, _modes=self.modes[self.mode_indices[i]:self.mode_indices[i + 1]]))

        return result


class Model_grid:
    """
    A grid of models.
    """

    def __init__(self):
        self.ndim = 0
        """
        Number of dimensions for the grid (excluding age), as based on the
        :py:data:`Model_grid.grid_params` variable
        """

        self.tracks = []
        """List of evolutionary tracks contained in the grid."""

        self.ndx = []
        """List containing track indices"""

        self.grid = None
        """Array containing the grid parameters for each evolutionary track (excluding age)."""

        self.tessellation = None
        """Object containing the tessellation of the grid used for interpolation."""

        self.grid_params = None
        """
        Set of parameters (excluding age) used to construct the grid and do interpolations.

        .. note::        
          For best interpolation results, these parameters should be comparable.
        """

        self.prefix = None
        """Root folder with grid of models (including final slash)."""

        self.postfix = ".freq"
        """Last part of the filenames which contain the model frequencies (default = ".freq")."""

        self.user_params = config.user_params
        """
        The set of user parameters involved in the grid.  This is to avoid having a different
        set of user parameters in `AIMS_configure.py`
        """

        self.distort_mat = None
        """
        Transformation matrix used to break the Cartesian character of the grid and reduce
        computation time
        """

    def read_model_list(self, filename):
        """
        Read list of models from a file and construct a grid.  If
        the mode format is Aldo, an entirely different strategy
        must be used.

        :param filename: name of the file with the list.
        :type filename: string
        """

        if (config.mode_format == "Aldo"):
            self.read_model_list_Aldo(filename)
        elif (config.mode_format == "BASTA"):
            self.read_model_list_BASTA(filename)
        else:
            self.read_model_list_standard(filename)

    def read_model_list_standard(self, filename):
        """
        Read list of models from a file and construct a grid.

        :param filename: name of the file with the list.  The first line
          of this file should contain a prefix which is typically the root
          folder of the grid of models.  This followed by a file with multiple
          columns.  The first 9 contain the following information for each model:

          1. the second part of the path.  When concatenated with the prefix
             on the first line, this should give the full path to the model.
          2. The stellar mass in :math:`\\mathrm{g}`
          3. The stellar radius in :math:`\\mathrm{cm}`
          4. The stellar luminosity in :math:`\\mathrm{g.cm^2.s^{-3}}`
          5. The metallicity
          6. The hydrogen content
          7. The stellar age in :math:`\\mathrm{Myrs}`
          8. The effective temperature in :math:`\\mathrm{K}`
          9. A dimensionless age parameter

          The following columns contain the parameters specified in the
          :py:data:`AIMS_configure.user_params` variable.

        :type filename: string
        """

        self.grid_params = config.grid_params

        # set the correct dimension:
        self.ndim = len(self.grid_params)

        # set prefix and postfix:
        with open(filename, "r") as listfile:
            line = listfile.readline().strip()
            columns = line.split()

            if len(columns) < 1:
                raise IOError("Erroneous first line in %s." % filename)

            self.prefix = columns[0]
            if (len(columns) > 1):
                self.postfix = columns[1]
            if (self.postfix == "None"):
                self.postfix = ""

        # read list with numpy:
        print("Reading list file")
        perm = [imass, iradius, iluminosity, iz0, ix0, iage, itemperature, iage_adim] \
               + [user_params_index[name_pair[0]] for name_pair in config.user_params]
        invperm = [1, ] * (len(perm) + 1)
        for i in range(len(perm)):
            invperm[perm[i]] = i + 1
        glb = np.loadtxt(filename, skiprows=1, usecols=invperm, dtype=gtype)
        names = np.loadtxt(filename, skiprows=1, usecols=(0,), dtype=filename_dtype)

        # sort list: this separates the evolutionary tracks (which still need
        #            to be partitionned)
        print("Sorting models")
        grid_list = [age_str, ] + list(self.grid_params)
        grid_temp = np.array([utilities.my_map(Model(glb[i]).string_to_param, \
                                               grid_list) for i in range(glb.shape[0])])
        ind = np.lexsort([grid_temp[:, i] for i in range(grid_temp.shape[1])])
        params_span = np.array([np.max(grid_temp[:, i]) - np.min(grid_temp[:, i]) \
                                for i in range(1, grid_temp.shape[1])])
        del grid_temp

        # sanity check:
        for i in range(len(self.grid_params)):
            span = params_span[i]
            if np.isnan(span):
                raise ValueError("The values of some models parameters are NaN.")
            if np.isinf(span):
                raise ValueError("The values of some models parameters are infinite.")
            if span == 0.0:
                raise ValueError("ERROR: parameter %s is constant in your grid" % self.grid_params[i] +
                                 "and cannot be used as a grid parameter.")

        # create tracks:
        print("Creating evolutionary tracks")
        nmodes = 0
        models_small_spectra = []
        for nmodels in range(glb.shape[0]):
            i = ind[nmodels]
            if isinstance(config.alpha_Fe_param, str):
                i_aFe = grid_list.index(config.alpha_Fe_param)
                aFe = glb[i][i_aFe]
            else:
                aFe = config.alpha_Fe_param
            aModel = Model(glb[i], _name=names[i], aFe=aFe)
            exceed_freqlim = aModel.read_file(self.prefix + names[i] + self.postfix)
            aModel.multiply_modes(1.0 / aModel.glb[ifreq_ref])  # make frequencies non-dimensional
            aModel.sort_modes()
            aModel.remove_duplicate_modes()
            if ((len(self.tracks) > 0) and (self.tracks[-1].matches(aModel))):
                self.tracks[-1].append(aModel)
            else:
                aTrack = Track(aModel, self.grid_params)
                self.tracks.append(aTrack)
            nmodes += len(aModel.modes)
            if (not exceed_freqlim):
                models_small_spectra.append(aModel.name)
            if (not config.batch):
                print("%d %d %d" % (len(self.tracks), nmodels + 1, nmodes))
                print('\033[2A')  # backup two line - might not work in all terminals
        print("%d %d %d" % (len(self.tracks), glb.shape[0], nmodes))
        del glb
        del names
        del ind

        # merge tracks which are too close:
        imerge = 0
        i = 1
        params_span *= eps
        while (i < len(self.tracks)):
            aModel = Model(self.tracks[i].glb[0], aFe=get_aFe(self.tracks[i]))
            for j in range(i):
                if (self.tracks[j].matches(aModel, params_tol=params_span)):
                    self.tracks[j].append_track(self.tracks[i])
                    del self.tracks[i]
                    imerge += 1
                    break
            else:
                i += 1
        print("I merged %d track(s)" % (imerge))

        # write list of models with spectra which are too small in a file:
        with open("models_small_spectra", "w") as output:
            for name in models_small_spectra:
                output.write(name + "\n")

        # sort tracks:
        for track in self.tracks:
            track.sort()

        # remove duplicate models:
        for track in self.tracks:
            if track.remove_duplicate_ages():
                print("WARNING: the track %s = %s" % (str(track.grid_params), str(track.params)))
                print("         has models with the same age.  Removing duplicate models.")

        # update list of indices:
        self.ndx = range(len(self.tracks))

        # need to create grid from scratch since tracks have been sorted.
        self.grid = np.asarray([track.params for track in self.tracks])

    def read_model_list_Aldo(self, filename):
        """
        Read list of models in Aldo format from a file and construct a grid.

        A track in Aldo format contains 3 parts:

          1. A three line header which is discarded.
          2. A set of four lines, for l=0 to 3, which describe the mode structure.  The
             first columns gives the start value of n, whereas the second column gives
             the number of modes for that particular l value.
          3. A table where each line corresponds to a model and the different columns
             correspond to the model ID, global quantities, and pulsation frequencies
             and inertias.

        :param filename: name of the file with the list of track files.
        :type filename: string
        """

        # sanity check:
        if (len(config.user_params) != 5):
            raise ValueError('user_params should have length 5 for Aldo format but got %i.' % len(config.user_params))

        for key in ['Xc', 'Zc', 'Xs', 'Zs', 'Mass_true']:
            if not key in user_params_index:
                raise KeyError("%s missing from user_params variable but needed for Aldo format." % key)

        self.grid_params = config.grid_params

        # set the correct dimension:
        self.ndim = len(self.grid_params)

        # set prefix and postfix:
        with open(filename, "r") as listfile:
            line = listfile.readline().strip()
            columns = line.split()

            if (len(columns) < 1):
                raise IOError("Erroneous first line in %s." % filename)

            self.prefix = columns[0]
            self.postfix = ""

            # read models and put them into evolutionary tracks:
            nmodels = 0
            nmodes = 0
            models_small_spectra = []
            norder = np.zeros((4, 2), dtype=ntype)
            for line in listfile:
                if line[0] == "#":
                    continue
                line = line.strip()
                if len(line) == 0:
                    continue

                # read track file:
                x0 = None
                z0 = None
                Mass0 = None
                with open(self.prefix + line, "r") as trackfile:
                    # skip header:
                    for i in range(3):
                        trackfile.readline()
                    # read mode structure
                    ntot = 0
                    for l in range(4):
                        columns = trackfile.readline().strip().split()
                        norder[l, 0] = int(columns[0])
                        norder[l, 1] = int(columns[1])
                        ntot += norder[l, 1]

                # read models with numpy:
                glbs = np.loadtxt(self.prefix + line, skiprows=7, usecols=range(1, 11 + 2 * ntot), dtype=gtype)
                names = np.loadtxt(self.prefix + line, skiprows=7, usecols=(0,), dtype=filename_dtype)

                # read models
                for i in range(glbs.shape[0]):
                    if (x0 is None):
                        x0 = glbs[i, 7]
                    if (z0 is None):
                        z0 = glbs[i, 6]
                    if (Mass0 is None):
                        Mass0 = glbs[i, 1] * constants.solar_mass
                    glb = np.empty((nglb,), dtype=gtype)
                    glb[imass] = Mass0
                    glb[iradius] = glbs[i, 4] * constants.solar_radius
                    glb[iluminosity] = constants.solar_luminosity * 10.0 ** glbs[i, 3]
                    glb[iz0] = z0
                    glb[ix0] = x0
                    glb[iage] = glbs[i, 0]
                    glb[itemperature] = glbs[i, 2]
                    glb[iage_adim] = glbs[i, 0]  # may be replace later on
                    glb[user_params_index["Mass_true"]] = glbs[i, 1]
                    glb[user_params_index["Xs"]] = glbs[i, 7]
                    glb[user_params_index["Zs"]] = glbs[i, 6]
                    glb[user_params_index["Xc"]] = glbs[i, 8]
                    glb[user_params_index["Zc"]] = 1.0 - glbs[i, 9] - glbs[i, 8]

                    aModel = Model(glb, _name=names[i], aFe=None)
                    exceed_freqlim = aModel.read_modes_Aldo(norder, glbs[i, 10:])
                    aModel.multiply_modes(1.0 / aModel.glb[ifreq_ref])  # make frequencies non-dimensional
                    aModel.sort_modes()
                    # We're assuming the models are grouped together in tracks, so
                    # we only need to check the previous track:
                    if len(self.tracks) > 0 and self.tracks[-1].matches(aModel):
                        self.tracks[-1].append(aModel)
                    else:
                        aTrack = Track(aModel, self.grid_params)
                        self.tracks.append(aTrack)
                    nmodels += 1
                    nmodes += len(aModel.modes)
                    if not exceed_freqlim:
                        models_small_spectra.append(aModel.name)
                    if not config.batch:
                        print("%d %d %d" % (len(self.tracks), nmodels, nmodes))
                        print('\033[2A')  # backup two line - might not work in all terminals

                del glbs
                del names

        print("%d %d %d" % (len(self.tracks), nmodels, nmodes))

        # find span for each grid parameter prior to merging tracks:
        grid_temp = np.asarray([track.params for track in self.tracks])
        params_span = np.array([np.max(grid_temp[:, i]) - np.min(grid_temp[:, i]) \
                               for i in range(grid_temp.shape[1])])
        del grid_temp

        # merge tracks which are too close:
        imerge = 0
        i = 1
        params_span *= eps
        while (i < len(self.tracks)):
            aModel = Model(self.tracks[i].glb[0], aFe=None)
            for j in range(i):
                if (self.tracks[j].matches(aModel, params_tol=params_span)):
                    self.tracks[j].append_track(self.tracks[i])
                    del self.tracks[i]
                    imerge += 1
                    break
            else:
                i += 1
        print("I merged %d track(s)" % (imerge))

        # write list of models with spectra which are too small in a file:
        with open("models_small_spectra", "w") as output:
            for name in models_small_spectra:
                output.write(name + "\n")

        # sort tracks:
        for track in self.tracks:
            track.sort()

        # sanity check:
        for track in self.tracks:
            if track.remove_duplicate_ages():
                print("WARNING: the track %s = %s" % (str(track.grid_params), str(track.params)))
                print("         has models with the same age.  Removing duplicate models.")

        # update list of indices:
        self.ndx = range(len(self.tracks))

        # need to create grid from scratch since tracks have been sorted.
        self.grid = np.asarray([track.params for track in self.tracks])

    def read_model_list_BASTA(self,filename):
        """
        Read list of models in BASTA format from a file and construct a grid.

        BASTA grids are contained HDF5 files with a tree-like structure.

        :param filename: name of the file with the list of track files.
        :type filename: string
        """

        # sanity check (to be completed if need be) ...

        qdict = {
            "FeH":"Fe_H",
            "FeHini":"Fe_H0",
            "MeH":"M_H",
            "MeHini":"M_H0",
            "alphaFe":"alpha_Fe",
            "alphaMLT":"alpha_MLT",
            "massini":"Mass0",
            "ove":"alpha_OV",
            "rho":"Rho",
            "rhocen":"RhoC",
            "xcen":"Xc",
            "xini":"X0",
            "xsur":"Xs",
            "ycen":"Yc",
            "yini":"Y0",
            "ysur":"Ys",
            "zcen":"Zc",
            "zini":"Z0",
            "zsur":"Zs",
        }

        qlatexdict = {
            "FeH":r"Iron content, $%s[Fe/H]%s$",
            "FeHini":r"Iron content, $%s[Fe/H]_0%s$",
            "LPhot":r"Luminosity, $%sL%s$",
            "Mbcz":r"Mass BCZ, $%sM_{\mathrm{BCZ}}%s$",
            "Mcore":r"Core mass, $%sM_{\mathrm{core}}%s$",
            "McoreX":r"H core mass, $%sM_{\mathrm{H,\,core}}%s$",
            "MeH":r"Metallicity, $%s[M/H]%s$",
            "MeHini":r"Metallicity, $%s[M/H]_0%s$",
            "Rbcz":r"Radius BCZ, $%sR_{\mathrm{BCZ}}%s$",
            "Rcore":r"Core radius, $%sR_{\mathrm{BCZ}}%s$",
            "RcoreX":r"H core radius, $%sR_{\mathrm{H,\,BCZ}}%s$",
            "Teff":r"Temperature, $%sT_{\mathrm{eff}}%s$",
            "ZAMSLPhot":r"ZAMS luminosity, $%sL_{\mathrm{ZAMS}}%s$",
            "ZAMSTeff":r"ZAMS temperature, $%sT_{\mathrm{eff,\,ZAMS}}%s$",
            "age":r"Age (Myrs)",
            "alphaFe":r"$%s\alpha_{\mathrm{Fe}}%s$",
            "alphaMLT":r"Mixing length parameter, $%s\alpha_{\mathrm{MLT}}%s$",
            "dnuscal":r"$%s\Delta\nu_{\mathrm{scal.}}%s$",
            "massfin":r"Mass, $%sM%s$",
            "massini":r"Initial mass, $%sM_0%s$",
            "numax":r"$%s\nu_{\mathrm{max}}%s$",
            "ove":r"Overshoot parameter, $%s\alpha_{\mathrm{ov.}}%s$",
            "radPhot":r"Photospheric radius, $%sR_{\mathrm{phot.}}%s$",
            "radTot":r"Total radius, $%sR_{\mathrm{tot.}}%s$",
            "rho":r"Mean density, $%s\rho%s$",
            "rhocen":r"Central density, $%s\rho_{\mathrm{C}}%s$",
            "tau0":r"Acoustic radius, $%s\tau%s$",
            "taubcz":r"BZC acoustic radius, $%s\tau_{\mathrm{BCZ}}%s$",
            "tauhe":r"He acoustic radius, $%s\tau_{\mathrm{He}}%s$",
            "xcen":r"Central hydrogen, $%sX_c%s$",
            "xini":r"Initial hydrogen, $%sX_0%s$",
            "xsur":r"Hydrogen content, $%sX_s%s$",
            "ycen":r"Central helium, $%sY_c%s$",
            "yini":r"Initial helium, $%sY_0%s$",
            "ysur":r"Helium content, $%sY_s%s$",
            "zcen":r"Central metallicity, $%sZ_c%s$",
            "zini":r"Initial metallicity, $%sZ_0%s$",
            "zsur":r"Metallicity, $%sZ_s%s$",
        }

        self.grid_params = config.grid_params

        # set the correct dimension:
        self.ndim = len(self.grid_params)

        # set prefix and postfix:
        self.prefix = "%s:"%(filename)
        self.postfix = ""

        # read grid
        h5_grid = h5py.File(filename)

        # read models and put them into evolutionary tracks:
        nmodels = 0
        nmodes  = 0

        # get names of quantities in tracks
        track0 = list(h5_grid["grid/tracks/"].keys())[0]
        quantities = list(h5_grid["grid/tracks/%s/"%(track0)].keys())
        quantities_del = ["name", "osc", "osckey", "volume_weight", "age", "massfin", \
                          "Teff", "zini", "xini", "radPhot", "LPhot"]
        for q in quantities_del:
            if (q in quantities):
                del quantities[quantities.index(q)]
        config.user_params = []
        for q in quantities:
            if (q in qdict):
                q1 = qdict[q]
            else:
                q1 = q
            if (q in qlatexdict):
                q2 = qlatexdict[q]
            else:
                q2 = q
            config.user_params.append((q1,q2))
        global nglb, nlin, ifreq_ref, iradius, iluminosity
        nglb        = 9 + len(quantities)
        nlin        = 6 + len(quantities)
        ifreq_ref   = 6 + len(quantities)
        iradius     = 7 + len(quantities)
        iluminosity = 8 + len(quantities)
        init_user_param_dict()
        self.user_params = config.user_params

        # read tracks
        for track in list(h5_grid["grid/tracks/"].keys()):

            # get models names:
            names = list(map(lambda s: "%s/%s"%(track,s.decode("utf-8")), \
                             h5_grid["grid/tracks/%s/name"%(track)]))
            nmodels_track = len(names)

            # get global quantities
            glb = np.zeros((nmodels_track,nglb),dtype=gtype)
            glb[:,iage]         = np.array(h5_grid["grid/tracks/%s/age"%(track)]) # (in Myrs ?)
            glb[:,imass]        = np.array(h5_grid["grid/tracks/%s/massfin"%(track)]) \
                                * constants.solar_mass # conversion to g
            glb[:,itemperature] = np.array(h5_grid["grid/tracks/%s/Teff"%(track)]) # K
            glb[:,iz0]          = np.array(h5_grid["grid/tracks/%s/zini"%(track)])
            glb[:,ix0]          = np.array(h5_grid["grid/tracks/%s/xini"%(track)])
            glb[:,iradius]      = np.array(h5_grid["grid/tracks/%s/radPhot"%(track)]) \
                                * constants.solar_radius # conversion to cm
            glb[:,iluminosity]  = np.array(h5_grid["grid/tracks/%s/LPhot"%(track)]) \
                                * constants.solar_luminosity # conversion to erg/s
            glb[:,ifreq_ref]    = 5e5*np.sqrt(constants.G*glb[:,imass]/glb[:,iradius]**3)/np.pi #muHz
            for q in quantities:
                if (q in qdict):
                    glb[:,user_params_index[qdict[q]]] = np.array(h5_grid["grid/tracks/%s/%s"%(track,q)])
                else:
                    glb[:,user_params_index[q]] = np.array(h5_grid["grid/tracks/%s/%s"%(track,q)])

            # initialise track
            aModel = Model(glb[0,:], _name=names[0], aFe=None)
            aTrack = Track(aModel,self.grid_params)
            aTrack.names = names
            aTrack.glb = glb

            # get mode frequencies and inertias
            nmodes_track = 0
            for n in range(nmodels_track):
                nmodes_track += len(h5_grid["grid/tracks/%s/osc"%(track)][n][0])
            aTrack.modes = np.empty((nmodes_track,),dtype=modetype)
            nmodes_track = 0
            aTrack.mode_indices = [nmodes_track,]
            for i in range(nmodels_track):
                n = len(h5_grid["grid/tracks/%s/osc"%(track)][i][0])
                aTrack.modes["n"][nmodes_track:nmodes_track+n]       = h5_grid["grid/tracks/%s/osckey"%(track)][i][1][:]
                aTrack.modes["l"][nmodes_track:nmodes_track+n]       = h5_grid["grid/tracks/%s/osckey"%(track)][i][0][:]
                aTrack.modes["freq"][nmodes_track:nmodes_track+n]    = h5_grid["grid/tracks/%s/osc"%(track)][i][0][:] \
                                                                     / glb[i,ifreq_ref]
                aTrack.modes["inertia"][nmodes_track:nmodes_track+n] = h5_grid["grid/tracks/%s/osc"%(track)][i][1][:]
                nmodes_track += n
                aTrack.mode_indices.append(nmodes_track)

            aTrack.nmodes = nmodes_track
            self.tracks.append(aTrack)
            nmodels += nmodels_track
            nmodes += nmodes_track
            del glb
            del names

            if (not config.batch):
                print("%d %d %d"%(len(self.tracks), nmodels, nmodes))
                print('\033[2A') # backup two line - might not work in all terminals

        h5_grid.close()
        print("%d %d %d"%(len(self.tracks), nmodels, nmodes))

        # find span for each grid parameter prior to merging tracks:
        grid_temp   = np.asarray([track.params for track in self.tracks])
        params_span = np.array([np.max(grid_temp[:,i])-np.min(grid_temp[:,i]) \
                               for i in range(grid_temp.shape[1])])
        del grid_temp

        # merge tracks which are too close:
        imerge = 0
        i = 1
        params_span *= eps
        while (i < len(self.tracks)):
            aModel = Model(self.tracks[i].glb[0], aFe=None)
            for j in range(i):
                if (self.tracks[j].matches(aModel,params_tol=params_span)):
                    self.tracks[j].append_track(self.tracks[i])
                    del self.tracks[i]
                    imerge += 1
                    break
            else:
                i+=1
        print("I merged %d track(s)"%(imerge))

        # sort tracks:
        for track in self.tracks: track.sort()

        # sanity check:
        for track in self.tracks:
            if track.remove_duplicate_ages():
                print("WARNING: the track %s = %s"%(str(track.grid_params),str(track.params)))
                print("         has models with the same age.  Removing duplicate models.")

        # update list of indices:
        self.ndx = range(len(self.tracks))

        # need to create grid from scratch since tracks have been sorted.
        self.grid = np.asarray([track.params for track in self.tracks])

    def replace_age_adim(self):
        """
        This replaces the dimensionless ages in the tracks according to the
        :py:data:`replace_age_adim` option chosen in ``AIMS_configure.py``.
        """

        if (config.replace_age_adim is None):
            return
        elif (config.replace_age_adim == "age"):
            for track in self.tracks:
                track.age_adim_to_age()
        elif (config.replace_age_adim == "scale_age"):
            for track in self.tracks:
                track.age_adim_to_scale_age()
        elif (config.replace_age_adim == "scale_Xc"):
            for track in self.tracks:
                track.age_adim_to_scale_Xc()
        else:
            raise ValueError("Unknown option %s for replace_scale_age in AIMS_configure.py." % config.replace_age_adim)

    def check_age_adim(self):
        """
        This checks that all of the tracks are sorted according to dimensionless age
        parameter.
        """

        for track in self.tracks:
            if not track.is_sorted_adim():
                raise ValueError("ERROR: track(s) not strictly sorted according to dimensionless age")

    def range(self, aParam):
        """
        Find range of values for the input parameter.

        :param aParam: name of the parameter for which to find the range
        :type aParam: str

        .. warning::
         The input parameter can only be one of the grid parameters or an
         age/mHe parameter.
        """

        if (aParam == age_str):
            param_min = self.tracks[0].age_lower
            param_max = self.tracks[0].age_upper

            for track in self.tracks:
                if (param_min > track.age_lower):
                    param_min = track.age_lower
                if (param_max < track.age_upper):
                    param_max = track.age_upper

        else:
            i = self.grid_params.index(aParam)
            param_min = self.tracks[0].params[i]
            param_max = self.tracks[0].params[i]

            for track in self.tracks:
                if (param_min > track.params[i]):
                    param_min = track.params[i]
                if (param_max < track.params[i]):
                    param_max = track.params[i]

        return [param_min, param_max]

    def tessellate(self):
        """Apply Delauny triangulation to obtain the grid tessellation."""

        if (self.distort_mat is None):
            self.tessellation = Delaunay(self.grid)
        else:
            self.tessellation = Delaunay(np.dot(self.grid, self.distort_mat))

    def plot_tessellation(self):
        """
        Plot the grid tessellation.

        .. warning::
          This only works for two-dimensional tessellations.
        """

        if (self.ndim != 2):
            print("Only able to plot the tessellation in two dimensions.")
            return

        # find bounds:
        xmin = np.nanmin(self.grid[:, 0])
        xmax = np.nanmax(self.grid[:, 0])
        ymin = np.nanmin(self.grid[:, 1])
        ymax = np.nanmax(self.grid[:, 1])
        dx = xmax - xmin
        dy = ymax - ymin
        xmin -= dx * 0.03
        xmax += dx * 0.03
        ymin -= dy * 0.05
        ymax += dy * 0.05

        # plt.semilogy(self.grid[:,0],self.grid[:,1],'o')
        plt.plot(self.grid[:, 0], self.grid[:, 1], 'o')
        plt.triplot(self.grid[:, 0], self.grid[:, 1], self.tessellation.simplices.copy())
        plt.xlim((xmin, xmax))
        plt.ylim((ymin, ymax))
        plt.xlabel(string_to_latex(self.grid_params[0]))
        plt.ylabel(string_to_latex(self.grid_params[1]))
        plt.savefig("tessellation.eps")

    def remove_tracks(self, nthreshold):
        """
        Removes stellar evolution tracks with fewer than nthreshold models.

        :param nthreshold: lower limit on number of models in a stellar
                           evolutionary track
        :type nthreshold: int

        :return: True if tracks have been removed and the grid needs to be
                 retessellated
        :rtype: boolean
        """

        # easy exit:
        if (nthreshold <= 0):
            return

        removedTracks = False
        for i in range(len(self.tracks) - 1, -1, -1):
            if (len(self.tracks[i].names) < nthreshold):
                del self.tracks[i]
                removedTracks = True

        # redo tessellation if need be
        if (removedTracks):
            print("I removed evolutionary tracks from the grid!")

            # update list of indices:
            self.ndx = range(len(self.tracks))

            # need to create grid from scratch since tracks have been sorted.
            self.grid = np.asarray([track.params for track in self.tracks])

        return removedTracks

    def distort_grid(self):
        """
        Define distortion matrix with which to distort grid to break its Cartesian
        character prior to tessellation.  This can cause find_simplex to run much
        faster.
        """

        self.distort_mat = np.dot(make_scale_matrix(self.grid), make_distort_matrix(self.ndim))
        print("Distortion matrix condition number: %e" % (np.linalg.cond(self.distort_mat)))

    def test_interpolation(self):
        """
        Test interpolation between different evolutionary tracks in a given grid.

        :return: The following four items are returned:

          - the interpolation errors
          - the first half of the partition (where the interpolation is tested)
          - the second half of the partition (used to carry out the interpolation)
          - the tessellation associated with the second half of the partition

        :rtype: np.array, list, list, tessellation object
        """

        ndx1, ndx2 = self.find_partition()
        if (self.distort_mat is None):
            tessellation = Delaunay(self.grid[ndx2, :])
        else:
            tessellation = Delaunay(np.dot(self.grid[ndx2, :], self.distort_mat))

        # initialisation
        results = []
        ndim = self.ndim + 1

        for j in ndx1:
            nmodels = len(self.tracks[j].names)
            aResult = np.empty((nmodels, ndim + nglb + 6), dtype=gtype)
            pt = list(self.tracks[j].params) + [0.0, ]

            for i in range(nmodels):
                aModel1 = Model(self.tracks[j].glb[i], aFe=None,
                                _modes=self.tracks[j].modes[self.tracks[j].mode_indices[i]: \
                                                            self.tracks[j].mode_indices[i + 1]])
                pt[-1] = aModel1.glb[iage]
                aModel2 = interpolate_model(self, pt, tessellation, ndx2)
                aResult[i, 0:ndim] = pt
                if (aModel2 is None):
                    aResult[i, ndim:ndim + nglb + 6] = np.nan
                else:
                    aResult[i, ndim:ndim + nglb + 6] = compare_models(aModel1, aModel2)

            results.append(aResult)

        return results, ndx1, ndx2, tessellation

    def find_partition(self):
        """
        Find a partition of the grid for use with :py:meth:`Model_grid.test_interpolation`

        :return: a random partition of [0 ... n-1] into two equal halves, where n is
                 the number of tracks in the grid
        :rtype: two lists of int

        """

        ndx = list(range(len(self.tracks)))
        random.shuffle(ndx)
        nn = len(self.tracks) // 2
        return ndx[:nn], ndx[nn:]

    def test_freq(self):
        """
        Test to see if frequencies in all of the models of the grid
        are in ascending order for each l value.

        :return: The following items are returned

          - the effective temperatures of the models with frequencies out of order
          - the luminosities of the models with frequencies out of order
          - the effective temperatures of the models with sorted frequencies
          - the luminosities of the models with sorted frequencies

        :rtype: four lists of floats
        """

        Teffs, Lums, Teffs_out, Lums_out = [], [], [], []
        for track in self.tracks:
            for i in range(len(track.names)):
                if (not track.freq_sorted(i)):
                    print(track.names[i])
                    Teffs_out.append(track.glb[i, itemperature])
                    Lums_out.append(math.log10(track.glb[i, iluminosity] / constants.solar_luminosity))
                else:
                    Teffs.append(track.glb[i, itemperature])
                    Lums.append(math.log10(track.glb[i, iluminosity] / constants.solar_luminosity))
        return Teffs_out, Lums_out, Teffs, Lums

    def find_epsilons(self, ltarget):
        """
        Find epsilon values in models from the grid

        :param ltarget: target l value for which epsilons are being obtained
        :type ltarget: int

        :return: the epsilon values
        :rtype: list of floats
        """

        epsilons = []
        for track in self.tracks:
            for i in range(len(track.names)):
                epsilon = track.find_epsilon(i, ltarget)
                if (epsilon != 0.0):
                    epsilons.append(epsilon)
        return epsilons


def init_user_param_dict():
    """
    Initialise the dictionaries which are related to user-defined parameters.  For
    a given parameter, these dictionaries provide the appropriate index for for
    the :py:data:`Model.glb` array as well as the appropriate latex name.
    """

    i = nlin - len(config.user_params)
    for (name, latex_name) in config.user_params:
        user_params_index[name] = i
        user_params_latex[name] = latex_name
        i += 1


def make_distort_matrix(d, theta=0.157):
    """
    Create a distortion matrix which can be used to make the grid more
    "tessellation-friendly", i.e. which leads to much shorter computation
    times for finding simplices.

    :param d: number of dimensions
    :type d: int

    :param theta: a small angle
    :type theta: float

    :return: a distortion matrix
    :rtype: 2D float array
    """
    cost = math.cos(theta)
    sint = math.sin(theta)
    mat = np.eye(d)
    aux = np.empty((d, d), dtype=float)
    for i in range(d):
        for j in range(i + 1, d):
            aux[:, :] = 0.0
            for k in range(d):
                if ((k == i) or (k == j)):
                    continue
                aux[k, k] = 1.0
            aux[i, i] = 1.0
            aux[j, j] = cost
            aux[j, i] = 0.0
            aux[i, j] = sint
            mat = np.dot(mat, aux)

    return mat


def make_scale_matrix(grid):
    """
    Create a distortion matrix which can be used to make the grid more
    "tessellation-friendly", i.e. which leads to much shorter computation
    times for finding simplices.

    :param grid: set of points used in the construction of the tessellation
    :type d: 2D float array

    :return: a distortion matrix
    :rtype: 2D float array
    """

    eps = 1e-3
    d = grid.shape[1]
    mat = np.eye(d)
    for i in range(d):
        x = np.asarray(list(set(grid[:, i])))
        x = np.sort(x)
        x = x[1:] - x[:-1]
        x = np.sort(x)
        ind = bisect_right(x, eps)
        mat[i, i] = 1.0 / x[ind]
        print("Scale factor %d: %e" % (i, mat[i, i]))
    return mat


def combine_models(model1, coef1, model2, coef2):
    """
    Do linear combination of this model with another.

    This method returns a new model which is the weighted sum
    of two models for the purposes of model interpolation.
    The classical parameters are combined in a self-consistent
    way as are the frequencies.

    :param model1: first model
    :param coef1:  weighting coefficient applied to first model
    :param model2: second model
    :param coef2:  weighting coefficient applied to second model

    :type model1: :py:class:`Model`
    :type coef1:  float
    :type model2: :py:class:`Model`
    :type coef2:  float

    :return: the combined model
    :rtype: :py:class:`Model`

    .. warning::
      One should avoid negative or zero coefficients as
      these could lead to undefined results.
    """

    # find global parameters (try to be self-consistent):

    # this first part is simply a linear combination:
    glb = np.empty((nglb,), dtype=gtype)
    glb[0:nlin] = coef1 * model1.glb[0:nlin] + coef2 * model2.glb[0:nlin]

    # this next part depends on previous results:
    glb[iradius] = (glb[imass] / (coef1 * model1.glb[imass] / model1.glb[iradius] ** 3
                                  + coef2 * model2.glb[imass] / model2.glb[iradius] ** 3)) ** (1.0 / 3.0)
    cnst1 = model1.glb[iluminosity] / (model1.glb[iradius] ** 2 * model1.glb[itemperature] ** 4)
    cnst2 = model2.glb[iluminosity] / (model2.glb[iradius] ** 2 * model2.glb[itemperature] ** 4)
    glb[iluminosity] = (coef1 * cnst1 + coef2 * cnst2) * glb[iradius] ** 2 * glb[itemperature] ** 4
    # glb[ifreq_ref] will be correctly defined when the Model() constructor is invoked

    # interpolate spectra:
    size3 = min(model1.modes.shape[0], model2.modes.shape[0])
    nvalues = np.empty((size3,), dtype=ntype)
    lvalues = np.empty((size3,), dtype=ltype)
    fvalues = np.empty((size3,), dtype=ftype)
    ivalues = np.empty((size3,), dtype=ftype)

    nvalues, lvalues, fvalues, ivalues, n3 = aims_fortran.combine_modes( \
        coef1, model1.modes['n'], model1.modes['l'], model1.modes['freq'], model1.modes['inertia'], \
        coef2, model2.modes['n'], model2.modes['l'], model2.modes['freq'], model2.modes['inertia'], \
        nvalues, lvalues, fvalues, ivalues)
    new_Model = Model(glb, _modes=list(zip(nvalues[0:n3], lvalues[0:n3], fvalues[0:n3], ivalues[0:n3])), aFe=None)
    new_Model.A_FeH = model1.A_FeH * coef1 + model2.A_FeH * coef2
    return new_Model


def compare_models(model1, model2):
    """
    Compare two models and find the largest frequency different for
    radial and non-radial modes.

    :param model1: first model
    :param model2: second model

    :type model1: :py:class:`Model`
    :type model2: :py:class:`Model`

    :return: a 1D array to be used in ``plot_test_interpolation.py``
      with the following measurements of the differences between the two models:

      - ``result[0]`` = maximum error on the radial modes
      - ``result[1]`` = RMS error on the radial modes
      - ``result[2]`` = RMS error on the radial modes near
        :math:`\\nu_{\\mathrm{max}}`
      - ``result[3]`` = maximum error on the non radial modes
      - ``result[4]`` = RMS error on the non radial modes
      - ``result[5]`` = RMS error on the non radial modes near
        :math:`\\nu_{\\mathrm{max}}`
      - ``result[6+[0:nglb]]`` = errors on the global parameters

    :rtype: np.array
    """

    if (config.interpolation_test_units == "microHz"):
        C1, C2 = model1.glb[ifreq_ref], model1.glb[ifreq_ref]
    elif (config.interpolation_test_units is None):
        C1 = C2 = 1.0
    else:
        raise ValueError("ERROR: Unrecognised units for interpolation_test_units in AIMS_configure.py")

    # initialisation:
    n_radial = 0
    n_radial_numax = 0
    n_non_radial = 0
    n_non_radial_numax = 0
    result = np.zeros((6 + nglb,), dtype=gtype)
    # define frequency interval around numax:
    numax = 0.5 * (C1 * model1.numax / model1.glb[ifreq_ref] \
                   + C2 * model2.numax / model2.glb[ifreq_ref])
    a = 0.8 * numax
    b = 1.2 * numax

    # compare frequency spectra:
    size1 = len(model1.modes)
    size2 = len(model2.modes)
    i1 = i2 = 0
    while ((i1 < size1) and (i2 < size2)):
        if (model1.modes['l'][i1] < model2.modes['l'][i2]):
            i1 += 1;
            continue
        if (model1.modes['l'][i1] > model2.modes['l'][i2]):
            i2 += 1;
            continue
        if (model1.modes['n'][i1] < model2.modes['n'][i2]):
            i1 += 1;
            continue
        if (model1.modes['n'][i1] > model2.modes['n'][i2]):
            i2 += 1;
            continue

        # now the two modes have the same n and l values:
        diff = abs(C1 * model1.modes['freq'][i1] - C2 * model2.modes['freq'][i2])
        avg_freq = (C1 * model1.modes['freq'][i1] + C2 * model2.modes['freq'][i2]) / 2.0
        if (model1.modes['l'][i1] == 0):
            if (result[0] < diff):
                result[0] = diff
            diff *= diff  # square diff
            result[1] += diff
            n_radial += 1
            # in python, this is called an interval comparison:
            if (a <= avg_freq <= b):
                result[2] += diff
                n_radial_numax += 1
        else:
            if (result[3] < diff):
                result[3] = diff
            diff *= diff  # square diff
            result[4] += diff
            n_non_radial += 1
            if (a <= avg_freq <= b):
                result[5] += diff
                n_non_radial_numax += 1
        i1 += 1
        i2 += 1

    # avoid divisions by zero:
    if (n_radial > 0):
        result[1] = math.sqrt(result[1] / float(n_radial))
    else:
        result[1] = np.nan

    if (n_radial_numax > 0):
        result[2] = math.sqrt(result[2] / float(n_radial_numax))
    else:
        result[2] = np.nan

    if (n_non_radial > 0):
        result[4] = math.sqrt(result[4] / float(n_non_radial))
    else:
        result[4] = np.nan

    if (n_non_radial_numax > 0):
        result[5] = math.sqrt(result[5] / float(n_non_radial_numax))
    else:
        result[5] = np.nan

    # absolute differences on global parameters:
    result[6:6 + nglb] = np.absolute(model1.glb - model2.glb)

    return result


def find_interpolation_coefficients(grid, pt, tessellation, ndx):
    """
    Find interpolation weights from the corresponding simplex.

    Linear interpolation weights are obtained with the simplex
    by finding the barycentric coordinates of the point given
    by ``pt``.

    :param grid: grid of models in which we're carrying out the
      interpolation
    :param pt: set of parameters used for finding the
      interpolation weights.  The first part contains the grid
      parameters (relevant to this interpolation), whereas
      the last element is the age (not used here).  If the
      provided set of parameters lies outside the grid, then
      ``None`` is returned instead of an interpolated model.
    :param tessellation: tessellation with which to carry out the
      interpolation.
    :param ndx: indices of the grid points associated with the
      tessellation

    :type grid: :py:class:`Model_grid`
    :type pt: array-like
    :type tessellation: Delaunay tessellation object
    :type ndx: list of int

    :return: lists of interpolation coefficients and tracks
    :rtype: list of floats, list of :py:class:`Track`
    """

    if (pt is None):
        return None, None
    if (tessellation is None):
        return None, None
    if (grid.distort_mat is not None):
        pt1 = np.dot(np.asarray(pt[0:-1], dtype=gtype), grid.distort_mat)
    else:
        pt1 = np.asarray(pt[0:-1], dtype=gtype)
    val = tessellation.find_simplex(pt1.reshape((1, grid.ndim)))[0]

    # see if point is outside tessellation
    if (val == -1):
        return None, None
    mat = tessellation.transform[val]

    # make sure the transformation matrix is defined:
    if (math.isnan(np.sum(mat))):
        return None, None

    b = mat[:grid.ndim].dot(pt1 - mat[grid.ndim])
    coefs = np.r_[b, 1.0 - b.sum()]
    ind = tessellation.simplices[val]

    # check to make sure you're not outside the grid:
    for coef in coefs:
        if (coef < -tol):
            return None, None

    # produce results, filtering out zero elements:
    coefs_out = []
    tracks = []
    for coef, i in zip(coefs, ind):
        # remove negative coefficients to avoid problems.
        if (coef > 0.0):
            coefs_out.append(coef)
            tracks.append(grid.tracks[ndx[i]])

    return coefs_out, tracks


def find_ages(coefs, tracks, age):
    """
    Find ages to which each track needs to be interpolated for a specified
    age.  The variable :py:data:`age_interpolation` in AIMS_configure.py
    decides between the following options:

    1. ``age_interpolation`` = ``age``: each track is simply interpolated
       to ``age``.
    2. ``age_interpolation`` = ``scale_age``: the age of each model along
       each evolutionary track, including the interpolated track, is
       linearly mapped onto the interval [0,1].  A dimensionless parameter
       ``eta`` is obtained by interpolating ``age`` onto the interval
       [0,1], using the linear transformation associated with the
       interpolated track.  Using the parameter eta, a corresponding age
       is obtained along each track.
    3. ``age_interpolation`` = ``age_adim``: each track is interpolated
       to the appropriate dimensionless age so that the interpolated
       physical age reproduces the input age.  This dimensionless age
       parameter is calculated via a fortran subroutine.

    .. figure:: ./figures/age_interpolation.*
      :figclass: align-center

      This diagram illustrates age interpolation for the two
      first options, namely ``age`` and ``scale_age``, and shows
      the advantages of selecting the latter.

    :param coefs: interpolation coefficients used to weight each track.
    :param tracks:  evolutionary tracks involved in the interpolation.
    :param age: target age for the output interpolated model.

    :type coefs: list of floats
    :type tracks:  list of :py:class:`Track`
    :type age: float

    :return: the relevant age for each track
    :rtype: list of floats

    .. note::
      - the interpolation coefficients should add up to 1.0
      - there should be as many tracks as interpolation coefficients.
    """

    assert (len(coefs) == len(tracks)), "Mismatch between len(coefs) and len(tracks)"

    if (config.age_interpolation == "age"):
        return [age] * len(coefs)

    elif (config.age_interpolation == "scale_age"):
        age_s = 0.0
        age_f = 0.0
        for coef, track in zip(coefs, tracks):
            age_s += coef * track.glb[0, iage]
            age_f += coef * track.glb[-1, iage]

        eta = (age - age_s) / (age_f - age_s)

        # check to see if the age lies within the interpolated track:
        if (eta < 0.0):
            return None
        if (eta > 1.0):
            return None

        ages = []
        for coef, track in zip(coefs, tracks):
            ages.append((1.0 - eta) * track.glb[0, iage] + eta * track.glb[-1, iage])

        return ages

    elif (config.age_interpolation == "age_adim"):

        coefs = np.array(coefs)
        ntracks = len(tracks)
        sze = np.array([track.glb.shape[0] for track in tracks], dtype=int)
        max_sze = np.max(sze)
        age_adim_array = np.zeros((max_sze, ntracks), dtype=gtype)
        age_array = np.zeros((max_sze, ntracks), dtype=gtype)
        weights = np.zeros((ntracks,), dtype=gtype)
        indices = np.zeros((ntracks, 3), dtype=int)
        age_adim = np.nan  # this is important to catch failures

        for i in range(ntracks):
            nmodels = sze[i]
            age_array[0:nmodels, i] = tracks[i].glb[:, iage]
            age_adim_array[0:nmodels, i] = tracks[i].glb[:, iage_adim]

        age_adim, indices, weights = aims_fortran.find_tau(age_adim_array, \
                                                           age_array, coefs, sze, age, age_adim, indices, weights)

        if (math.isnan(age_adim)):
            return None

        ages = np.zeros((ntracks,), dtype=gtype)
        for i in range(ntracks):
            ages[i] = weights[i] * age_array[indices[i, 0] - 1, i] \
                      + (1.0 - weights[i]) * age_array[indices[i, 0], i]

        return ages

    else:

        raise ValueError("Unrecognised age_interpolation option %s." % config.age_interpolation)


def interpolate_model(grid, pt, tessellation, ndx):
    """
    Interpolate model in grid using provided parameters.

    The interpolation is carried out in two steps.  First, linear
    interpolation according to age is carried out on each node of
    the simplex containing the set of parameters.  This interpolation
    is done using the :py:class:`Track.interpolate_model` method.
    Then, linear interpolation is carried out within the simplex.
    This achieved by finding the barycentric coordinates of the
    model (i.e. the weights), before combining the age-interpolated
    models form the nodes using the :py:class:`combine_models` method.
    In this manner, the weights are only calculated once, thereby
    increasing computational efficiency.

    :param grid: grid of models in which we're carrying out the
      interpolation
    :param pt: set of parameters used for the interpolation.
      The first part contains the grid parameters, whereas
      the last element is the age.  If the provided set
      of parameters lies outside the grid, then ``None``
      is returned instead of an interpolated model.
    :param tessellation: tessellation with which to carry out the
      interpolation.
    :param ndx: indices of the grid points associated with the
      tessellation

    :type grid: :py:class:`Model_grid`
    :type pt: array-like
    :type tessellation: Delaunay tessellation object
    :type ndx: list of int

    :return: the interpolated model
    :rtype: :py:class:`Model`
    """

    # find simplex interpolation coefficients
    coefs, tracks = find_interpolation_coefficients(grid, pt, tessellation, ndx)
    if (coefs is None):
        return None

    # find ages:
    ages = find_ages(coefs, tracks, pt[-1])
    if (ages is None):
        return None

    n = len(tracks)

    # treat the case where there is only 1 model:
    if (n == 1):
        if (abs(coefs[0] - 1.0) > eps):
            print("WARNING: erroneous interpolation coefficient: " + str(coefs[0]))
        return tracks[0].interpolate_model(ages[0])

    # treat the case where there are at least 2 models:
    aModel1 = tracks[0].interpolate_model(ages[0])
    if (aModel1 is None):
        return None
    aModel2 = tracks[1].interpolate_model(ages[1])
    if (aModel2 is None):
        return None
    aModel1 = combine_models(aModel1, coefs[0], aModel2, coefs[1])
    for i in range(2, n):
        aModel2 = tracks[i].interpolate_model(ages[i])
        if (aModel2 is None):
            return None
        aModel1 = combine_models(aModel1, 1.0, aModel2, coefs[i])

    return aModel1


def find_combination(grid, pt):
    """
    Find linear combination of models which corresponds to interpolating
    the model based on the provided parameters.

    The interpolation is carried out using the same procedure as in
    :py:func:`interpolate_model`.

    :param grid: grid of models in which we're carrying out the
      interpolation
    :param pt: set of parameters used for the interpolation.
      The first part contains the grid parameters, whereas
      the last element is the age.  If the provided set
      of parameters lies outside the grid, then ``None``
      is returned instead of an interpolated model.

    :type grid: :py:class:`Model_grid`
    :type pt: array-like

    :return: pairs of coefficients and model names
    :rtype: tuple of (float,string)
    """

    # find simplex interpolation coefficients
    coefs, tracks = find_interpolation_coefficients(grid, pt, grid.tessellation, grid.ndx)
    if (coefs is None):
        return None

    # find ages:
    ages = find_ages(coefs, tracks, pt[-1])
    if (ages is None):
        return None

    n = len(tracks)

    # combine multiple models:
    results = ()
    for coef, track, age in zip(coefs, tracks, ages):
        if (coef < 0.0):  # make sure we're not outside the grid
            return None
        result = track.find_combination(age, coef)
        if (result is None):
            return None
        results += result
    return results
