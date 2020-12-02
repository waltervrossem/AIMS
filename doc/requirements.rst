Requirements
============

The following python packages are needed for AIMS:

  * `dill <https://pypi.python.org/pypi/dill/>`_
  * `emcee <http://dan.iel.fm/emcee/current/>`_
  * `corner <https://github.com/dfm/corner.py>`_

    - note: this used to be called triangle in previous releases

  * `numpy <http://www.numpy.org/>`_
  * `matplotlib <http://matplotlib.org/>`_
  * `multiprocessing <https://docs.python.org/2/library/multiprocessing.html>`_
  * `f2py <https://github.com/pearu/f2py/wiki>`_

    - note: this is usually included with numpy

  * `tqdm <https://pypi.org/project/tqdm/>`_

.. warning::
  As of version 3, emcee no longer supports parallel tempering.
  Accordingly, AIMS is only compatible with older versions.  One can specify
  a version of emcee during installation with the following command (we recommend
  using version 2.2.1)::

    sudo pip install emcee==2.2.1
