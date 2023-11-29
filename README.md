# PyRacmacs
A python interface to Racmacs using rpy2.

Almost all of the functionalities in Racmacs (and some more) are also available in PyRacmacs. Most properties and functions in Racmacs which would be in camelCase 
is in snake_case in PyRacmacs but otherwise available in a similar fashion. See examples for extensive tutorials on usage of PyRacmacs.

## Installation notes:
If you want to get this running in a conda environment but have R packages installed in your base, the only way for the environment to see the R packages installed 
in the base seems to be installing rpy2 via pip and setting the environment variables correctly before doing so. See:

https://stackoverflow.com/questions/68936589/how-to-select-r-installation-when-using-rpy2-on-conda

After making sure you have r-base installed in your enviroment, run the following:

```
export LDFLAGS="-Wl,-rpath,/usr/lib/R/lib" (change path to where base R lib is)
pip install rpy2 --force-reinstall --compile --no-binary rpy2
```

After the installation is complete, you can check by running python -m rpy2.situation. If somewhere in the output it says

```
Looking for R's HOME:
    Environment variable R_HOME: None
    Calling `R RHOME`: /usr/lib/R
```

you are good to go. Note that re-installing other R related packages such as r-base after this might modify this (export R_HOME = /usr/lib/R does not remedy this), in which case
you need to reinstall rpy2.

## Warning and error notes:
Some of the warning and error messages returned will directly be from R. An example is when there are disconnected points in a map, you will get the following when you run PyRacmacs optimize on such a map:

```
RRuntimeError: Error: Map contains disconnected points (points that are not connected through any path of detectable titers so cannot be coordinated relative to each other). To optimize anyway, rerun with 'options = list(ignore_disconnected = TRUE)'.
```

In such a case you should check for the corresponding way of supplying options in PyRacmacs. In PyRacmacs supplying options is done via dictionaries (as named lists don't exist in python) and the corresponding fix would be
supplying to the PyRacmacs optimize function

```
options={"ignore_disconnected":True}
```
