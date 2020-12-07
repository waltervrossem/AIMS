#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Warrick Ball
#    2020 December 7
#    Birmingham, UK
#

"""Start of a unit test framework for AIMS functionality."""

import pytest
import functions

################################################################################
################################################################################
################################################################################

delta = 1e-6

def test_identity():
    assert functions.identity([1.0]) == 1.0

def test_identity_gradient():
    numerical = 0.5*(functions.identity([1.0+delta])-functions.identity([1.0-delta]))/delta
    analytic = functions.identity_gradient([1.0])[0]
    assert pytest.approx(analytic, numerical)

def test_ratio():
    assert functions.ratio([1.0, 2.0]) == 0.5

def test_ratio_gradient():
    numerical = [0.5*(functions.ratio([1.0+delta, 2.0])-functions.ratio([1.0-delta, 2.0]))/delta,
                 0.5*(functions.ratio([1.0, 2.0+delta])-functions.ratio([1.0, 2.0-delta]))/delta]
    analytic = functions.ratio_gradient([1.0, 2.0])
                         
    assert pytest.approx(numerical[0], analytic[0])
    assert pytest.approx(numerical[1], analytic[1])

def test_norm():
    assert functions.norm([1.0, 2.0]) == 5.0**0.5

def test_norm_gradient():
    numerical = [0.5*(functions.norm([1.0+delta, 2.0])-functions.norm([1.0-delta, 2.0]))/delta,
                 0.5*(functions.norm([1.0, 2.0+delta])-functions.norm([1.0, 2.0-delta]))/delta]
    analytic = functions.norm_gradient([1.0, 2.0])
                         
    assert pytest.approx(numerical[0], analytic[0])
    assert pytest.approx(numerical[1], analytic[1])
