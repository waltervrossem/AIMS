#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Warrick Ball
#    2020 December 7
#    Birmingham, UK
#

"""Start of a unit test framework for AIMS functionality."""

import pytest
import utilities

################################################################################
################################################################################
################################################################################

def test_to_float():
    pytest.approx(utilities.to_float('1.0D-0'), 1.0)

def test_is_number():
    assert utilities.is_number('1.0d+0')
    assert not utilities.is_number('zzz')

def test_trim():
    assert utilities.trim('data # comment # comment') == 'data '
    
