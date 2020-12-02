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

"""
A module which contains various reference functions which are useful for frequency combination functions.
"""

import math

def identity(vec):
    """
    Apply the identity function to the input value (x).

    We recall that the identify function is expressed as follows:
    :math:`f(x) = x`

    :param vec: the input vector
    :type vec: float array like

    .. warning::
       No verifications are carried out within the function so the user
       must make sure the input vector only has one component, as
       implicitely assumed.
    """

    return vec[0]

def identity_gradient(vec):
    """
    Apply the gradient of the identity function to the input value (x).

    We recall that the gradient of the identify function is expressed as follows:
    :math:`\\vec{\\Delta} f(x) = 1`

    :param vec: the input vector
    :type vec: float array like

    .. warning::
       No verifications are carried out within the function so the user
       must make sure the input vector only has one component, as
       implicitely assumed.
    """

    return (1.0,)

def ratio(vec):
    """
    Return the ratio for input values (x,y).

    We recall that the ratio function is expressed as follows:
    :math:`f(x,y) = \\frac{x}{y}`

    :param vec: the input vector
    :type vec: float array like

    .. warning::
       No verifications are carried out within the function so the user
       must make sure the input vector only has two components, as
       implicitely assumed.
    """

    return vec[0]/vec[1]

def ratio_gradient(vec):
    """
    Return the gradient of the ratio function for input values (x,y).

    We recall that gradient of the ratio function is expressed as follows:
    :math:`\\vec{\\Delta} f(x,y) = \\left(\\frac{1}{y},-\\frac{x}{y^2}\\right)`

    :param vec: the input vector
    :type vec: float array like

    .. warning::
       No verifications are carried out within the function so the user
       must make sure the input vector only has two components, as
       implicitely assumed.
    """

    return (1.0/vec[1],-vec[0]/(vec[1]*vec[1]))

def norm(vec):
    """
    Return the norm of an input vector.

    We recall that the norm is expressed as follows:
    :math:`\\|(\\vec{v})\\| = \\sqrt{\\sum_i v_i^2}`

    :param vec: the input vector
    :type vec: float array like
    """

    return math.sqrt(sum([e*e for e in vec]))

def norm_gradient(vec):
    """
    Return the gradient of the norm function for an input vector.

    We recall that the gradient of the norm function takes on the following
    expression:
    :math:`\\vec{\\Delta}\\|(\\vec{v})\\| = \\frac{\\vec{v}}{\\|\\vec{v}\\|}`

    :param vec: the input vector
    :type vec: float array like
    """

    my_norm = math.sqrt(sum([e*e for e in vec]))
    return [e/my_norm for e in vec] 
