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
    d = 5
    mat = model.make_distort_matrix(d, theta=0.0)
    for i in range(d):
        for j in range(d):
            if i == j:
                assert mat[i,j] == 1.0
            else:
                assert mat[i,j] == 0.0

def test_find_ages():
    test = model.Model_grid()
    test.read_model_list('tests/data/test.aimslist')

    t = 1000.
    age_interpolation = model.config.age_interpolation
    for option in ['age', 'scale_age', 'age_adim']:
        ages = model.find_ages([1.0], test.tracks[:1], t)
        assert ages[0] == pytest.approx(t)

    model.config.age_interpolation = age_interpolation

def test_Model():
    test = model.Model_grid()
    test.read_model_list('tests/data/test.aimslist')
    m = test.tracks[0].interpolate_model(1000.0)
    
    # m.read_file_agsm('tests/data/modelS.agsm') # never passed
    
    m.write_file_simple('tests/data/tmp.simple')
    m.read_file_CLES('tests/data/tmp.simple') # doubles up all modes
    m.sort_modes()
    assert m.remove_duplicate_modes() # docs say opposite of what happens

    m.append_modes(m.modes[-1])
    assert not m.remove_duplicate_modes()
    assert m.remove_duplicate_modes()

    assert m.get_age() == pytest.approx(1000.0)

    assert m.get_freq()[0] == m.modes['freq'][0]
    assert m.find_mode(m.modes['n'][0], m.modes['l'][0]) == m.modes['freq'][0]
    nmin, nmax, lmin, lmax = m.find_mode_range()
    assert nmin <= nmax
    assert lmin == 0
    assert lmax == 0

    assert 0.5 < m.numax/m.cutoff < 0.7 # vaguely reasonable values

    assert m.freq_sorted()

    m.print_me()

def test_Track():
    test = model.Model_grid()
    test.read_model_list('tests/data/test.aimslist')

    track = test.tracks[0]

    assert track.is_sorted()
    assert track.freq_sorted(0)
    assert not track.duplicate_ages()
    assert not track.remove_duplicate_ages()
    m = track.interpolate_model(1000.0)

def test_Model_grid():
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
    assert 0 <= t[1] <= 1e12
    
    assert pytest.approx(test.range('Z'), [0.018, 0.022])

    test.test_interpolation()

    Teffs_out, Lums_out, Teffs, Lums = test.test_freq()
    # everything should be in order
    assert len(Teffs_out) == 0
    assert len(Lums_out) == 0

    epsilons = test.find_epsilons(0)
    assert 0 < epsilons[0] < 2

    assert test.remove_tracks(20) # should be all of them
    
