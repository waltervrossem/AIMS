Requirements
============

The following python packages are needed for AIMS:

  * `dill <https://pypi.org/project/dill/>`_
  * `emcee <https://emcee.readthedocs.io/en/stable/>`_ 
  * `ptemcee <https://github.com/waltervrossem/ptemcee>`_

    - this is only needed if you decide to do parallel tempering
    - note: This is a different version than the one available on PyPI as this version contains a bugfix for newer versions of python. It can be installed with
    ``pip install ptemcee@git+https://github.com/waltervrossem/ptemcee@f6e91e7``

  * `corner <https://corner.readthedocs.io/en/latest/>`_

    - note: this used to be called triangle in previous releases

  * `scipy <https://scipy.org/>`_
  * `numpy <https://numpy.org/>`_
  * `f2py <https://numpy.org/doc/stable/f2py/>`_
  
    - this is usually included with numpy

  * `matplotlib <https://matplotlib.org/>`_
  * `multiprocessing <https://docs.python.org/3/library/multiprocessing.html>`_
  
    - this is already part of the standard python library

  * `tqdm <https://pypi.org/project/tqdm/>`_
  * `lxml <https://lxml.de/>`_

For convenience, a `requirements.txt` file has been included.  This allows
the user to install the needed python packages via the command::

    pip install -r requirements.txt

.. note::
  This version of SPInS is compatible both with versions 2 and 3 of emcee,
  an MCMC package written in python (see `Foreman-Mackey et al., 2013, PASP
  125, 306 <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_).

  However, as of version 3, emcee no longer supports parallel tempering.
  Therefore, SPInS now uses the separate package ptemcee for parallel
  tempering. This also has the added benefit of enabling the use of
  dynamic temperatures which may speed up convergence in some cases
  (see `Vousden, Farr, and Mandel, 2016, MNRAS 455, 1919
  <https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.1919V/abstract>`_)
