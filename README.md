# PyRacmacs
A python interface to Racmacs using rpy2.

Almost all of the functionalities in Racmacs (and some more) are also available in PyRacmacs. Most properties and functions in Racmacs which would be in camelCase 
is in snake_case in PyRacmacs but otherwise available in a similar fashion. See examples for extensive tutorials on usage of PyRacmacs.

## Installation notes:
If you want to get this running in a conda environment but have R packages installed in your base, the only way for the environment to see the R packages installed 
in the base seems to be installing rpy2 via pip and setting the environment variables correctly before doing so. See:

https://stackoverflow.com/questions/68936589/how-to-select-r-installation-when-using-rpy2-on-conda

A minimal working example is:

1- Create an enviroment containing pip
```
conda create -n env-name pip 
conda activate env-name
```

2- Check that which R points to the base R and which R points to env bin
```
which R
usr/bin/R
which pip
path_to_env-name/bin/pip
```

3- export the necessary LDFLAG and install rpy2 via pip from your enviroment:
```
export LDFLAGS="-Wl,-rpath,/usr/lib/R/lib"
pip install rpy2 --force-reinstall --compile --no-binary rpy2
```

You can check if the installation was correct by running python -m rpy2.situation.
If somewhere in the output it says

```
Looking for R's HOME:
    Environment variable R_HOME: None
    Calling `R RHOME`: /usr/lib/R
```

you are good to go. 

Note that installing other R related packages such as r-base prior to this will result in R RHOME likely pointing to enviroment's bin directory
and I am not sure how to create such an environment, apart from trying to install r-base afterwards rpy2.

If you want to create an environment with more packages changing the first step to 

```
conda create -n env-name r-base other-packages
```

seems to work fine too.


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
