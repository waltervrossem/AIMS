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
    
def test_sparse_print():
    tmpfile = 'tmp.sparse_print'
    utilities.sparse_print(tmpfile, utilities.np.eye(4))
    with open(tmpfile, 'r') as f:
        for i, line in enumerate(f.readlines()):
            pytest.approx(utilities.to_float(line.split()[0]), 1.0)
            assert(line.split()[1] == '(%i,' % i)
            assert(line.split()[2] == '%i)' % i)
