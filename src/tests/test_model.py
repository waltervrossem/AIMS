#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Warrick Ball
#    2020 December 7
#    Birmingham, UK
#

"""Start of a unit test framework for AIMS functionality."""

import pytest
import model

################################################################################
################################################################################
################################################################################

def test_init_user_param_dict():
    assert model.init_user_param_dict() == None

def test_make_distort_matrix():
    """
    Test the creation of a trivial distortion matrix.
    """
    d = 5
    mat = model.make_distort_matrix(d, theta=0.0)
    for i in range(d):
        for j in range(d):
            if i == j:
                assert mat[i,j] == 1.0
            else:
                assert mat[i,j] == 0.0

def test_make_scale_matrix():
    """
    Test the creation of a trivial scale matrix.
    """
    grid = model.np.eye(5)
    mat = model.make_scale_matrix(grid)
    assert mat == pytest.approx(grid)

def test_find_ages():
    """
    Test the find_ages function.
    """
    test = model.Model_grid()
    test.read_model_list('tests/data/test.aimslist')

    t = 1000.
    age_interpolation = model.config.age_interpolation

    for option in ['age', 'scale_age', 'age_adim']:
        model.config.age_interpolation = option
        ages = model.find_ages([1.0], test.tracks[:1], t)
        assert ages[0] == pytest.approx(t)

    model.config.age_interpolation = age_interpolation

def test_Model():
    """
    Test functions in the Model class.
    """
    test = model.Model_grid()
    test.read_model_list('tests/data/test.aimslist')
    m = test.tracks[0].interpolate_model(1000.0)

    # choose options through model.config to cover more code
    # so save original option to restore after tests
    mode_format = model.config.mode_format
    model.config.mode_format = 'agsm'
    m.read_file('tests/data/test.agsm')
    assert all(m.modes['n'] == [19,20,21])
    assert all(m.modes['l'] == 1)
    assert m.modes['freq'] == pytest.approx(
        [2830.96281472,2967.04152439,3102.93893379])
    assert m.modes['inertia'] == pytest.approx(
        [3.91089245e-10,3.18748640e-10,2.64390334e-10])

    m.write_file_simple('tests/data/tmp.simple')
    model.config.mode_format = 'simple'
    m.read_file('tests/data/tmp.simple')
    model.config.mode_format = 'CLES'
    m.read_file('tests/data/tmp.simple')
    model.config.mode_format = 'CLES_Mod'
    m.read_file_CLES_Mod('tests/data/test.cles_mod')
    model.config.mode_format = 'PLATO'
    m.read_file_PLATO('tests/data/test.plato')
    assert not m.freq_sorted() # by design of the plato test data
    m.sort_modes()
    assert m.freq_sorted()

    model.config.mode_format = mode_format

    m.append_modes(m.modes[-1])            # duplicate a mode
    assert not m.remove_duplicate_modes()  # removing it should return False
    assert m.remove_duplicate_modes()      # nothing to remove should return True

    assert m.get_age() == pytest.approx(1000.0) # age should be what we interpolated to

    for option in [None, 'Kjeldsen2008', 'Kjeldsen2008_scaling',
                   'Kjeldsen2008_2', 'Ball2014', 'Ball2014_2',
                   'Sonoi2015', 'Sonoi2015_scaling', 'Sonoi2015_2']:
        assert all(m.get_freq(a=[0.0, 0.0, 0.0], surface_option=option)
                   == m.modes['freq'])

    assert m.find_mode(m.modes['n'][0], m.modes['l'][0]) == m.modes['freq'][0]
    assert model.np.isnan(m.find_mode(99, 99))
    nmin, nmax, lmin, lmax = m.find_mode_range()
    assert nmin <= nmax
    assert lmin == 0
    assert lmax == 2

    assert 0.5 < m.numax/m.cutoff < 0.7 # vaguely reasonable values

    assert 0 < m.find_epsilon(0) < 2
    assert 0 == m.find_epsilon(99)

    c = model.combine_models(m, 1.0, m, 0.0)
    comparison = model.compare_models(m, m)
    assert comparison[0] == pytest.approx(0.0)
    assert comparison[1] == pytest.approx(0.0)
    assert comparison[3] == pytest.approx(0.0)
    assert comparison[4] == pytest.approx(0.0)
    assert comparison[6:10] == pytest.approx([0.0, 0.0, 0.0, 0.0])

    # test log, ln and exp transformations using a few different
    # variables and check that L ~ R²Teff⁴
    assert m.string_to_param('log_Mass') == model.math.log10(m.string_to_param('Mass'))
    assert m.string_to_param('Radius') == model.math.exp(m.string_to_param('ln_Radius'))
    assert m.string_to_param('Luminosity') == model.math.log(m.string_to_param('exp_Luminosity'))
    assert m.string_to_param('Luminosity') == pytest.approx(
        m.string_to_param('Radius')**2*(m.string_to_param('Teff')/5777.)**4, rel=1e-3) # uncertainty because of Teff_sun

    assert sum([m.string_to_param(k) for k in 'XYZ']) == pytest.approx(1.0) # abundances must sum to 1

    # surface abundances don't change in test grid
    assert m.FeH == pytest.approx(m.FeH0)
    assert m.MH == pytest.approx(m.MH0)
    assert m.zsx_s == pytest.approx(m.zsx_0)

    m.modes = model.np.empty(0, dtype=model.modetype)
    assert m.find_mode_range() == (-1, -1, -1, -1)

    m.print_me()
    del(m)

def test_Track():
    """
    Test functions in the Track class.
    """
    test = model.Model_grid()
    test.read_model_list('tests/data/test.aimslist')

    track = test.tracks[0]

    assert track.is_sorted()
    assert track.freq_sorted(0)
    # reverse the modes in one model so they definitely aren't sorted
    track.modes[track.mode_indices[0]:track.mode_indices[1]] = \
        track.modes[track.mode_indices[0]:track.mode_indices[1]][::-1]
    assert not track.freq_sorted(0)
    assert not track.duplicate_ages()
    assert not track.remove_duplicate_ages()
    m = track.interpolate_model(1000.0)

    c = track.find_combination(1000.0, 1.0)
    assert c[0][0] + c[1][0] == pytest.approx(1.0)

    ages, freqs = track.find_modes(10, 0)
    assert ages[0] <= ages[-1]

    ages_dim, freqs_dim = track.find_modes_dim(10, 0)
    assert ages[0]/ages[1] == pytest.approx(ages_dim[0]/ages_dim[1])

    nmin, nmax, lmin, lmax = track.find_mode_range()
    assert nmin <= nmax
    assert lmin == 0
    assert lmax == 2

    assert track.age_upper-track.age_lower == pytest.approx(track.age_range)

    # if we append the track to itself, there should be duplicates
    track.append_track(track)
    assert track.remove_duplicate_ages()

def test_Model_grid():
    """
    Test functions in the Model_grid class.
    """
    test = model.Model_grid()
    test.read_model_list('tests/data/test.aimslist')

    replace_age_adim = model.config.replace_age_adim
    # 'scale_Xc' is also an option but Xc isn't in the grid
    for option in [None, 'age', 'scale_age']:
        model.config.replace_age_adim = option
        test.replace_age_adim()

    model.config.replace_age_adim = replace_age_adim

    test.check_age_adim()
    t = test.range('Age')
    assert 0 <= t[0]
    assert 0 <= t[1] <= 20e9

    assert pytest.approx(test.range('Z'), [0.018, 0.022])

    test.test_interpolation()

    Teffs_out, Lums_out, Teffs, Lums = test.test_freq()
    # everything should be in order
    assert len(Teffs_out) == 0
    assert len(Lums_out) == 0

    epsilons = test.find_epsilons(0)
    assert 0 < epsilons[0] < 2
    assert test.find_epsilons(99) == []

    assert test.remove_tracks(20) # should be all of them
