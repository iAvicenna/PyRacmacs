# PyRacmacs
A python interface to Racmacs using rpy2.

Almost all of the functionalities in Racmacs (and some more) are also available in PyRacmacs. Most properties and functions in Racmacs which would be in camelCase 
is in snake_case in PyRacmacs but otherwise available in a similar fashion. See examples for extensive tutorials on usage of PyRacmacs.

Installation notes:
If you want to get this running in a conda environment but have R packages installed in your base, the only way for the environment to see the R packages installed 
in the base seems to be installing rpy2 via pip and setting the environment variables correctly before doing so. See:

https://stackoverflow.com/questions/68936589/how-to-select-r-installation-when-using-rpy2-on-conda

1- export LDFLAGS="-Wl,-rpath,/usr/lib/R/lib" (change path to where base R lib is)

2- pip install rpy2 --force-reinstall --compile --no-binary rpy2

After the installation is complete, you can check by running python -m rpy2.situation. If somewhere in the output it says

Looking for R's HOME:
    Environment variable R_HOME: None
    Calling `R RHOME`: /usr/lib/R

Note that installing other R related packages such as r-base might modify this (export R_HOME = /usr/lib/R does not remedy this), in which case
you need to reinstall rpy2.
