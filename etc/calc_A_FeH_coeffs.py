#!/usr/bin/env python
# coding: utf-8
# $Id: AIMS.py
# Author: Walter E. van Rossem <walter.vanrossem@unibo.it>
# Copyright (C) Walter E. van Rossem and contributors
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
Script to generate A_FeH conversion factors using https://github.com/aarondotter/initial_xa_calculator
"""

import os
import tqdm
import numpy as np
from scipy import optimize as op

# Assumes AIMS and initial_xa_calculator are in the same directory
xa_calc = '../../initial_xa_calculator/initial_xa_calculator'
net_name = 'pp_and_cno_extras.net'
feh = 0.0
DYDZ = 2

alpha_Fe_list = np.linspace(-4, 1, 101)

metallicities = []
for aFe in tqdm.tqdm(alpha_Fe_list):
    out = os.popen(f'/home/walter/Github/initial_xa_calculator/initial_xa_calculator {net_name} {feh} {aFe} {DYDZ} && rm input_initial_xa.data input_XYZ').read()
    out = [float(_.split('=')[1]) for _ in out.split('\n')[1:3]]  # ([M/H], [Fe/H])
    metallicities.append(out)
metallicities = np.array(metallicities)

def func(aFe, a, b):
    return np.log10(a * 10**aFe + b)

popt, pcov = op.curve_fit(func, alpha_Fe_list, metallicities[:,0])
print(f'A_FeH = {popt[0]} * 10**[a/Fe] + {popt[1]}')
