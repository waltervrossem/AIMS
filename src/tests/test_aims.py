#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Steven Hale
#    2020 May 8
#    Birmingham, UK
#

"""Start of a unit test framework for AIMS functionality."""

import pytest
import AIMS

################################################################################
################################################################################
################################################################################

def test_configuration():
    """
    Run the check_configuration() function to make sure the
    configuration variables are acceptable.
    """

    # Raises AssertionError on failure.
    AIMS.check_configuration()

################################################################################

def test_Distribution():
   """
   Test creating an instance of the Distribution class.
   """

   test = AIMS.Distribution("Gaussian", 1)

def test_Prior_list():
   """
   Test creating an instance of the Prior_list class.
   """

   test = AIMS.Prior_list()

def test_Mode():
   """
   Test creating an instance of the Mode class.
   """

   test = AIMS.Mode(1, 1, 1, 1)

def test_Combination():
   """
   Test creating an instance of the Combination class.
   """

   test = AIMS.Combination()

def test_Likelihood():
   """
   Test creating an instance of the Likelihood class.
   """

   test = AIMS.Likelihood()

def test_Probability():
   """
   Test creating an instance of the Probability class.
   """

   test = AIMS.Probability(1, 1)

################################################################################
