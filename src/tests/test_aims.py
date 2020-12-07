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

   # distributions are currently only defined up to an additive constant
   # so we only test differences in case the constant changes

   test = AIMS.Distribution("Gaussian", [0.0, 1.0])
   assert test.mean == 0.0
   assert test.error_bar == 1.0
   assert test(0.0) == test(1.0) + 0.5

   # third parameter is how many sigma to truncate
   test = AIMS.Distribution("Truncated_gaussian", [0.0, 1.0, 5.0])
   assert test.mean == 0.0
   assert test.error_bar == 1.0
   assert test(0.0) == test(1.0) + 0.5

   test = AIMS.Distribution("Uniform", [0.0, 1.0]) # unit uniform
   assert test.mean == 0.5
   assert test.error_bar == 0.5
   assert test(0.0) == test(1.0)

   test = AIMS.Distribution("Uninformative", [])
   assert AIMS.np.isnan(test.mean)
   assert AIMS.np.isnan(test.error_bar)
   assert test(0.0) == test(1.0)

   test = AIMS.Distribution("Above", [0.0])
   assert AIMS.np.isnan(test.mean)
   assert AIMS.np.isnan(test.error_bar)
   assert test(0.0) == test(1.0)

   test = AIMS.Distribution("Below", [0.0])
   assert AIMS.np.isnan(test.mean)
   assert AIMS.np.isnan(test.error_bar)
   assert test(0.0) == test(-1.0)

def test_Prior_list():
   """
   Test creating an instance of the Prior_list class.
   """

   test = AIMS.Prior_list()
   assert test([]) == 0.0
   test.add_prior(AIMS.Distribution("Gaussian", [0.0, 1.0]))
   assert test([0.0]) == test([1.0]) + 0.5

def test_Mode():
   """
   Test creating an instance of the Mode class.
   """

   mode1 = AIMS.Mode(1, 1, 1, 1)
   mode2 = AIMS.Mode(1, 1, 1, 1)
   assert mode1.match(mode2)

def test_Combination():
   """
   Test creating an instance of the Combination class.
   """
   test = AIMS.Combination()
   test.add_coeff(0, 0.5)
   test.add_coeff(1, 0.5)
   assert test.find_values([0.0, 1.0]) == 0.5

def test_Combination_function():
   """
   Test creating an instance of the Combination_function class.
   """
   # numerator
   c_num = AIMS.Combination()
   c_num.add_coeff(0, 1.0)
   c_num.add_coeff(1, 0.0)

   # denominator
   c_den = AIMS.Combination()
   c_den.add_coeff(0, 0.0)
   c_den.add_coeff(1, 1.0)

   test = AIMS.Combination_function('test',
                                    lambda x: x[0]/x[1],
                                    lambda x: (1.0/x[1], -x[0]/x[1]**2))
   test.add_combination(c_num)
   test.add_combination(c_den)
   assert test.find_values([1.0, 2.0]) == 0.5

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

def test_write_grid():
    """
    Tests the creation of a binary grid from plain text data.
    """

    AIMS.write_binary_data('tests/data/test.aimslist',
                           'tests/data/test.aimsgrid')

def test_interpolation_test():
    """
    Tests the interpolation test.
    """

    test = AIMS.interpolation_tests('tests/data/interpolation_test')

def test_grid_functions():
    grid = AIMS.load_binary_data('tests/data/test.aimsgrid')
    assert grid.range('Z') == [0.018, 0.022]
    assert grid.tracks[0].is_sorted()
    assert grid.tracks[0].is_sorted_adim()
    assert grid.tracks[0].freq_sorted(0)

def test_fit_data():
    """
    Tests fitting the models to data.  Based on AIMS script.
    """
    config = AIMS.config
    AIMS.check_configuration()
    like = AIMS.Likelihood()
    like.read_constraints('tests/data/test.aimsobs', factor=1.0)
    like.guess_dnu(with_n=True)
    AIMS.utilities.my_map(like.add_seismic_constraint,config.seismic_constraints)
    like.find_covariance()
    like.create_combination_arrays()
    like.find_weights()

    grid = AIMS.load_binary_data(config.binary_grid)
    grid_params_MCMC = grid.grid_params + (AIMS.model.age_str,)
    grid_params_MCMC_with_surf = grid_params_MCMC \
        + AIMS.model.get_surface_parameter_names(config.surface_option)
    nsurf            = len(AIMS.model.get_surface_parameter_names(config.surface_option))
    ndims            = len(grid_params_MCMC_with_surf)
    nwalkers         = 4*ndims # something small for fast testing

    AIMS.ndims = ndims
    AIMS.nsurf = nsurf

    priors = AIMS.Prior_list()

    for param_name in grid_params_MCMC_with_surf:
        priors.add_prior(AIMS.Distribution("Uninformative",[]))

    prob = AIMS.Probability(priors,like)

    pool = None
    my_map = AIMS.utilities.my_map

    p0 = AIMS.np.random.uniform(low =[0.999, 0.0199, 4500.0],
                                high=[1.001, 0.0201, 5000.0],
                                size=(nwalkers, ndims))

    sampler = AIMS.emcee.EnsembleSampler(nwalkers, ndims, prob)
    AIMS.grid = grid
    pos, prob, state = sampler.run_mcmc(p0, 200)
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos, 200)

################################################################################
