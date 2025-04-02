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
   Test various functions in the Distribution class.
   """

   # distributions are currently only defined up to an additive constant
   # so we only test differences in case the constant changes

   test = AIMS.Distribution("Gaussian", [0.0, 1.0])
   assert test.mean == 0.0
   assert test.error_bar == 1.0
   assert test(0.0) == test(1.0) + 0.5
   test.realisation() # could be anywhere
   assert test.nparams == 2

   test.re_centre(1.0)
   assert test.mean == 1.0

   test.re_normalise(2.0)
   assert test.error_bar == 2.0

   # third parameter is how many sigma to truncate
   test = AIMS.Distribution("Truncated_gaussian", [0.0, 1.0, 5.0])
   assert test.mean == 0.0
   assert test.error_bar == 1.0
   assert test(0.0) == test(1.0) + 0.5
   assert test(-10.0) == test(10.0) # same value outside support
   assert -5.0 <= test.realisation() <= 5.0
   assert test.nparams == 3

   test.re_centre(1.0)
   assert test.mean == 1.0

   test.re_normalise(2.0)
   assert test.error_bar == 2.0

   test = AIMS.Distribution("Uniform", [0.0, 1.0]) # unit uniform
   assert test.mean == 0.5
   assert test.error_bar == 0.5
   assert test(0.0) == test(1.0)
   assert test(-1.0) == test(2.0) # same value outside support
   assert 0.0 <= test.realisation() <= 1.0
   assert test.nparams == 2

   test.re_centre(1.0)
   assert test.mean == 1.0

   test.re_normalise(2.0)
   assert test.error_bar == 2.0

   test = AIMS.Distribution("Uninformative", [])
   assert AIMS.np.isnan(test.mean)
   assert AIMS.np.isnan(test.error_bar)
   assert test(0.0) == test(1.0)
   assert test.nparams == 0

   with pytest.raises(ValueError):
       test.realisation()

   with pytest.raises(ValueError):
       test.re_normalise(1.0)

   test = AIMS.Distribution("Above", [0.0])
   assert AIMS.np.isnan(test.mean)
   assert AIMS.np.isnan(test.error_bar)
   assert test(0.0) == test(1.0)
   assert test(-1.0) == test(-2.0) # same value outside support
   assert test.nparams == 1

   with pytest.raises(ValueError):
       test.realisation()

   with pytest.raises(ValueError):
       test.re_normalise(1.0)

   test = AIMS.Distribution("Below", [0.0])
   assert AIMS.np.isnan(test.mean)
   assert AIMS.np.isnan(test.error_bar)
   assert test(0.0) == test(-1.0)
   assert test(1.0) == test(2.0) # same value outside support
   assert test.nparams == 1

   with pytest.raises(ValueError):
       test.realisation()

   with pytest.raises(ValueError):
       test.re_normalise(1.0)

   # IMF1 is a power law between param[0] and param[1] with index param[2]
   # IMF2 is a two-part power law with index param [3] between param[0] and param[1] and index param[4] between param[1] and param[2]
   IMF1 = AIMS.Distribution("IMF1", [0.5, 1.5, -1.0])
   IMF2 = AIMS.Distribution("IMF2", [0.5, 1.0, 1.5, -1.0, -1.0])

   assert IMF1(0.0) == IMF1(2.0)
   assert 0.5 <= IMF1.realisation() <= 1.5
   assert IMF1.nparams == 3
   with pytest.raises(ValueError):
       IMF1.re_centre(0.0)

   with pytest.raises(ValueError):
       IMF1.re_normalise(1.0)

   assert IMF2(0.0) == IMF2(2.0)
   assert 0.5 <= IMF2.realisation() <= 1.5
   assert IMF2.nparams == 5
   with pytest.raises(ValueError):
       IMF2.re_centre(0.0)

   with pytest.raises(ValueError):
       IMF2.re_normalise(1.0)

   assert IMF1(0.8)-IMF1(1.2) == pytest.approx(IMF2(0.8)-IMF2(1.2))
   assert IMF1.mean == pytest.approx(IMF2.mean)
   assert IMF1.error_bar == pytest.approx(IMF2.error_bar)

   test.print_me()
   test.to_string()

   test = AIMS.Distribution("Undefined", [1, 2, 3])
   assert test.nparams == 0

   with pytest.raises(ValueError):
       test(0.0)

   with pytest.raises(ValueError):
       test.mean

   with pytest.raises(ValueError):
       test.error_bar

   with pytest.raises(ValueError):
       test.realisation()

   with pytest.raises(ValueError):
       test.re_centre(0.0)

   with pytest.raises(ValueError):
       test.re_normalise(1.0)

def test_Prior_list():
   """
   Test creating an instance of the Prior_list class.
   """

   test = AIMS.Prior_list()
   assert test([]) == 0.0
   test.add_prior(AIMS.Distribution("Uniform", [0.0, 1.0]))
   assert test([0.0]) == test([1.0])
   assert 0.0 <= test.realisation() <= 1.0
   assert (0.0 <= test.realisation(size=2)).all()
   assert (1.0 >= test.realisation(size=2)).all()
   assert (0.0 <= test.realisation(size=(2,2))).all()
   assert (1.0 >= test.realisation(size=(2,2))).all()

def test_Mode():
   """
   Test the Mode class.
   """

   mode1 = AIMS.Mode(1, 1, 1, 1)
   mode2 = AIMS.Mode(1, 1, 1, 1)
   assert mode1.match(mode2)

   mode2 = AIMS.Mode(1, 2, 3, 4)
   assert not mode1.match(mode2)

   mode1.print_me()

def test_Combination():
   """
   Tests the Combination class by creating a Combination
   that represents (x[0]+x[1])/2.
   """
   test = AIMS.Combination()
   test.add_coeff(0, 0.5)
   test.add_coeff(1, 0.5)
   assert test.find_values([0.0, 1.0]) == 0.5

   test.print_me()

def test_Combination_function():
   """
   Test the Combination_function class.
   """
   # numerator: given [x[0], x[1]], returns x[0]
   c_num = AIMS.Combination()
   c_num.add_coeff(0, 1.0)
   c_num.add_coeff(1, 0.0)

   # denominator: given [x[0], x[1]], returns x[1]
   c_den = AIMS.Combination()
   c_den.add_coeff(0, 0.0)
   c_den.add_coeff(1, 1.0)

   # combination is x[0]/x[1]
   test = AIMS.Combination_function('test',
                                    lambda x: x[0]/x[1],
                                    lambda x: (1.0/x[1], -x[0]/x[1]**2))
   test.add_combination(c_num)
   test.add_combination(c_den)
   assert test.find_values([1.0, 2.0]) == 0.5

   test.print_me()

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

    AIMS.model.init_user_param_dict()
    AIMS.write_binary_data('tests/data/test.aimslist',
                           'tests/data/test.aimsgrid')

def test_interpolation_test():
    """
    Tests the interpolation test.
    """

    test = AIMS.interpolation_tests('tests/data/interpolation_test')

def test_grid_functions():
    AIMS.grid = AIMS.load_binary_data('tests/data/test.aimsgrid')
    assert AIMS.grid.range('Z') == pytest.approx([0.018, 0.022])
    assert AIMS.grid.tracks[0].is_sorted()
    assert AIMS.grid.tracks[0].is_sorted_adim()
    assert AIMS.grid.tracks[0].freq_sorted(0)

    AIMS.write_list_file('tests/data/tmp.aimslist')
    AIMS.write_SPInS_file_cgs('tests/data/tmp.spins')
    AIMS.write_SPInS_file_solar('tests/data/tmp.spins')

def test_fit_data():
    """
    Tests fitting the models to data.  Based on AIMS script.
    Tests a variety of other things along the way.
    """
    # we change some settings in this test, so need to save them
    config = AIMS.config
    AIMS.check_configuration()
    AIMS.my_map = AIMS.utilities.my_map

    AIMS.like = AIMS.Likelihood()
    AIMS.like.clear_seismic_constraints()
    AIMS.like.read_constraints('tests/data/test.aimsobs', factor=1.0)
    AIMS.my_map(AIMS.like.add_seismic_constraint,
                ["r02","r01","r10", "avg_dnu","avg_dnu0", "skip",
                 "dnu", "dnu0", "d2nu", "nu_min0","nu_min1","nu_min2"])
    AIMS.like.clear_seismic_constraints()

    AIMS.my_map(AIMS.like.add_seismic_constraint, ["nu"])
    AIMS.like.find_covariance()
    AIMS.like.clear_seismic_constraints()

    AIMS.like.add_dnu_constraint_matrix("dnu", l_targets=[])
    AIMS.like.add_dnu_constraint_matrix("dnu")
    AIMS.like.add_dnu_constraint("dnu", l_targets=[])
    AIMS.like.clear_seismic_constraints()

    AIMS.like.guess_dnu(with_n=True)
    AIMS.my_map(AIMS.like.add_seismic_constraint, config.seismic_constraints)
    AIMS.like.find_covariance()
    AIMS.like.create_combination_arrays()

    weight_option = config.weight_option

    for option in [None, "Absolute", "Relative", "Reduced",
                   "Reduced_bis"]:
        config.weight_option = option
        AIMS.like.find_weights()

    config.weight_option = "aoj178nsdfja12"
    with pytest.raises(ValueError):
        AIMS.like.find_weights()

    config.weight_option = weight_option
    AIMS.like.find_weights()

    AIMS.config.user_params = () # causes to be overwritten from grid
    grid = AIMS.load_binary_data(config.binary_grid)
    grid_params_MCMC = grid.grid_params + (AIMS.model.age_str,)
    grid_params_MCMC_with_surf = grid_params_MCMC \
        + AIMS.model.get_surface_parameter_names(config.surface_option)
    nsurf            = len(AIMS.model.get_surface_parameter_names(config.surface_option))
    ndims            = len(grid_params_MCMC_with_surf)
    nwalkers         = 4*ndims # something small for fast testing

    AIMS.grid_params_MCMC = grid_params_MCMC
    AIMS.grid_params_MCMC_with_surf = grid_params_MCMC_with_surf
    AIMS.ndims = ndims
    AIMS.nsurf = nsurf
    AIMS.config.nwalkers = nwalkers

    AIMS.priors = AIMS.Prior_list()

    for param_name in grid_params_MCMC_with_surf:
        AIMS.priors.add_prior(AIMS.Distribution("Uninformative",[]))

    AIMS.prob = AIMS.Probability(AIMS.priors, AIMS.like)

    # config.tight_ball_range is a dict of initial distributions
    # key is parameter name; value[1] is width
    # we steadily build up to 5x observed constraints on "Mass" and "Z"
    # so that MCMC will narrow distribution
    # testing code branches on the way to final initial distribution
    AIMS.config.tight_ball_range = {}
    p0 = AIMS.init_walkers()
    AIMS.config.tight_ball_range["Mass"] = ("Gaussian", [0.0, 0.005])
    p0 = AIMS.init_walkers()
    AIMS.config.tight_ball_range["Z"] = ("Gaussian", [0.0, 0.0005])
    p0 = AIMS.init_walkers()

    assert AIMS.prob.is_outside(-p0[0])
    assert not AIMS.prob.is_outside(p0[0])

    # the best model should be the same one from which we took the
    # obsevational data (tests/data/0.020_1.00_7)
    AIMS.find_best_model()
    assert AIMS.best_grid_result > -0.1 # should be very good
    # assert AIMS.best_grid_params[0] == pytest.approx(1.0) # fails because of difference in Msun
    assert AIMS.best_grid_params[1] == pytest.approx(0.02)      # Z
    assert AIMS.best_grid_params[2] == pytest.approx(4772.9476) # Age
    assert AIMS.best_grid_params[5] == pytest.approx(5761.1637) # Teff

    mode_map, nmissing = AIMS.like.find_map(AIMS.best_grid_model, True)
    mode_map, nmissing = AIMS.like.find_map(AIMS.best_grid_model, False)

    nsurf = AIMS.nsurf
    surface_option = config.surface_option
    for option in ["Kjeldsen2008", "Kjeldsen2008_2",
                   "Kjeldsen2008_scaling", "Ball2014", "Ball2014_2",
                   "Sonoi2015", "Sonoi2015_scaling", "Sonoi2015_2"]:
        config.surface_option = option
        AIMS.nsurf = len(AIMS.model.get_surface_parameter_names(option))
        AIMS.like.get_optimal_surface_amplitudes(AIMS.best_grid_model,
                                                 mode_map)
    AIMS.nsurf = nsurf
    config.surface_option = surface_option

    AIMS.output_folder = 'tests/data'
    AIMS.write_model(AIMS.best_grid_model, AIMS.best_grid_params,
                     AIMS.best_grid_result, "tmp",
                     extended=True)

    # we hit a few more code branches by using a surface correction
    surface_option = AIMS.config.surface_option
    AIMS.config.surface_option = "Sonoi2015"
    AIMS.config.output_osm = 'tests/data'
    AIMS.write_osm_frequencies('tmp.osm_freq', AIMS.best_grid_model)
    AIMS.write_osm_don('tmp.osm_don', AIMS.best_grid_model)
    AIMS.write_osm_xml('tmp.osm_xml', AIMS.best_grid_params,
                       AIMS.best_grid_model)
    AIMS.config.surface_option = surface_option

    sampler = AIMS.emcee.EnsembleSampler(nwalkers, ndims, AIMS.prob)
    AIMS.grid = grid
    pos, prob, state = sampler.run_mcmc(p0, 200)
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos, 200)

    # quantitative contraints are 10x observations in data/test.aimsobs
    assert sampler.get_chain(flat=True)[:,0].mean() == pytest.approx(1.0, abs=0.01)
    assert sampler.get_chain(flat=True)[:,1].mean() == pytest.approx(0.02, abs=0.001)

    # all this really tests is that the output has narrowed the
    # distribution compared to the initial uniform distribution p0,
    # which was deliberately broad
    assert sampler.get_chain(flat=True)[:,0].std() < p0[:,0].std()
    assert sampler.get_chain(flat=True)[:,1].std() < p0[:,1].std()

    AIMS.write_samples('tests/data/tmp.samples', ['Mass'], sampler.get_chain(flat=True)[:,:1])
    AIMS.write_statistics('tests/data/tmp.stats', ['Mass'], sampler.get_chain(flat=True)[:,:1])
    AIMS.write_percentiles('tests/data/tmp.percentiles', ['Mass'], sampler.get_chain(flat=True)[:,:1])
    AIMS.config.tight_ball = False
    AIMS.write_readme('tests/data/tmp.readme', 0.0)
    AIMS.write_combinations('tests/data/tmp.combinations', sampler.get_chain(flat=True))

    # garbage data in the correct format
    LEGACY_keys = ['Radius', 'Mass', 'log_g', 'Rho', 'Age', 'Teff',
                   'Fe_H', 'Luminosity', 'X', 'Y', 'Xc']
    AIMS.write_LEGACY_summary('tests/data/tmp.legacy', '0', LEGACY_keys,
                              sampler.get_chain(flat=True)[:,[0]*len(LEGACY_keys)])

    # restore settings
    AIMS.config = config

################################################################################
