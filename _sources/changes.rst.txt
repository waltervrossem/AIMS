List of changes
===============

+--------------+---------------------------------------------------------------------------+
| **Version**  | **Changes**                                                               |
+--------------+---------------------------------------------------------------------------+
| 2.2.0        | * fixed compilation error for compatibility with python>=3.12 and numpy>=2|
|              | * add option to plot a subset of params in corner plot                    |
|              | * add option to make writing samples files optional                       |
|              | * docs use version in ``AIMS.py``                                         |
|              | * empty output folder before generating output                            |
|              | * use ``venv`` in CI and get CI base ``requirements.txt`` from            |
|              |   src/requirements.txt                                                    |
|              | * Add option to add 'non-constraints' which are only used during plotting |
|              | * use ``AIMS_configure.py`` from current working directory and a copy of  |
|              |   it and the input file are saved in the output folder                    |
|              | * fixed various deprecation warnings                                      |
|              | * use a version of ptemcee with a fix for indexing error in newer versions|
|              |   of numpy                                                                |
|              | * fixed bug when reading gyre files which skipped the first line of data  |
|              | * use [M/H] calculation for [Fe/H]                                        |
|              | * also plot nu%dnu + dnu in echelle diagrams                              |
|              | * replace remaining sys.exit with errors                                  |
+--------------+---------------------------------------------------------------------------+
| 2.1.0        | * added compatibility with emcee3 and ptemcee                             |
|              | * expanded unit test coverage                                             |
|              | * added ability to run AIMS without seismic constraints                   |
|              | * various bugfixes                                                        |
+--------------+---------------------------------------------------------------------------+
| 2.0.0        | * modified storage for evolutionary tracks thus saving a lot of memory    |
|              | * a dimensionless age parameter is introduced for the purposes of         |
|              |   interpolation between tracks (using the same approach as in SPInS)      |
|              | * new and more flexible implementation of frequency combinations          |
|              | * added a batch mode (without status bar)                                 |
+--------------+---------------------------------------------------------------------------+
| 1.3.0        | * can run both in python2.x and python3.x                                 |
|              | * included default distributions for priors and tight ball ranges         |
|              | * can tolerate erroneous values of user_params when fitting observations  |
|              | * added status bar for iterations (thanks to B. Rendle)                   |
|              | * added plots with evolution of walker percentiles                        |
+--------------+---------------------------------------------------------------------------+
| 1.2.0        | * removed extrapolation beyond grid limits                                |
|              | * various subprograms rewritten in Fortran (thus accelerating the code)   |
+--------------+---------------------------------------------------------------------------+
| 1.1.0        | * added extrapolation beyond grid limits                                  |
+--------------+---------------------------------------------------------------------------+
| 1.0.0        | * initial version                                                         |
+--------------+---------------------------------------------------------------------------+
