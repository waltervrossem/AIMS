List of changes
===============

+--------------+---------------------------------------------------------------------------+
| **Version**  | **Changes**                                                               |
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
