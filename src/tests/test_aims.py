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
