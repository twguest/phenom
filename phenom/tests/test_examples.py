import numpy as np

from phenom.spectra import flat_spectra


def test_one_plus_one_is_two():
    "Check that one and one are indeed two."
    assert 1 + 1 == 2


def test_spectra():
    assert np.sum(flat_spectra()) == 50
