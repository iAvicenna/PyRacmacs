# PyRacmacs
A python interface to Racmacs using rpy2.

Almost all of the functionalities in Racmacs (and some more) are also available in PyRacmacs. Most properties and functions in Racmacs which would be in camelCase 
is in snake_case in PyRacmacs but otherwise available in a similar fashion. See examples for extensive tutorials on usage of PyRacmacs.

Installation notes:
If you want to get this running in a conda environment but have R packages installed in your base, the only way for the environment to see the R packages installed 
in the base seems to be installing rpy2 via pip. Before doing it run "which pip" to make sure it points to the envrionment you want to install this in. Follow
the instructions here:

https://stackoverflow.com/questions/51486081/install-and-use-rpy2-using-conda-so-that-it-uses-default-r-installation-in-us

If it is still not working you can try to diagnose if the installation was correct (such as what Rpath rpy2 uses) by running "python -m rpy2.situation"

Also see here:
https://unix.stackexchange.com/questions/637298/arch-linux-arm-error-while-loading-shared-libraries-libblas-so-3-cannot-open
