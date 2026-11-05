Installation
============

To avoid issues with other python installations, it is recommended to use a virtual
environment and install the required packages there, activate the environment, and
then compile it::

    python -m venv .venv --prompt AIMS
    source .venv/bin/activate
    pip install AIMS@git+https://github.com/waltervrossem/AIMS

then to use first activate the python environment and then run ::

    source .venv/bin/activate
    AIMS observation_file

To install the development version and environment use::

    pip install AIMS@git+https://github.com/waltervrossem/AIMS@dev --group dev

