---
title: Galileo100 @ CINECA
updated: 2022-09-27 15:12:48Z
created: 2022-09-20 07:03:54Z
latitude: 45.46542190
longitude: 9.18592430
altitude: 0.0000
---

# python3 environment

create your own python environment in your $HOME
```
module load python/3.8.12--gcc--10.2.0
python3 -m venv myenv
source $HOME/myenv/bin/activate
```

modify the activate script by adding on top the following lines 
```
module load python/3.8.12--gcc--10.2.0
module load py-setuptools/58.2.0--gcc--10.2.0
module load py-numpy/1.21.3--gcc--10.2.0
```

re-source the environment 
```
deactivate
source $HOME/myenv/bin/activate
```

and install ipython for convenience 
```
pip3 install ipython
``` 

# WarpX

## Build

## Run

## Post-process

create your own python virtual environment
pip3 install ipython

# Smilei
## Build
download the source code in your `$HOME`
```
git clone https://github.com/SmileiPIC/Smilei.git smilei
```

move to the right directory
```
cd smilei
``` 


## Run

## Post-process

# EPOCH
## Build
clone the repository
```
git clone --recursive https://github.com/Warwick-Plasma/epoch.git
```

create a `epoch.profile` file in your `$HOME` with the following lines
```
module purge
module load intel/oneapi-2022--binary
module load intelmpi/oneapi-2022--binary
module load python/3.8.12--intel--2021.4.0
```

then move to the target directory, e.g.
```
cd epoch/epoch2d
```

and compile by typing
```
source $HOME/epoch.profile
make -j8 COMPILER=intel
```

the executable should be generated in `epoch/epoch2d/bin`

# Run
create a directory in $CINECA_SCRATCH
cd $CINECA_SCRATCH
mkdir epoch_tests

copy an input file in this directory 

# Post-process EPOCH data
source your python3 environment 
move to any epoch source directory (e.g. `epoch/epoch2d`) and type
```
make sdfutils
```
a messagge similar to this should appear
```
writing list of installed files to 'pybuild/files.txt'
```
you have to copy all the files listed in `/epoch/SDF/utilities/pybuild/files.txt` in the correct directory of your python environment `$HOME/myenv/lib/python3.8/site-packages/`
this should work
```
export MYDIR=$HOME/myenv/lib/python3.8/site-packages/
cp /g100/home/userexternal/aforment/.local/lib/python3.8/site-packages/sdf_legacy.py $MYDIR
cp -r /g100/home/userexternal/aforment/.local/lib/python3.8/site-packages/sdf_helper $MYDIR/sdf_helper
cp /g100/home/userexternal/aforment/.local/lib/python3.8/site-packages/__pycache__/sdf_legacy.cpython-38.pyc $MYDIR/__pycache__
cp /g100/home/userexternal/aforment/.local/lib/python3.8/site-packages/sdf.cpython-38-x86_64-linux-gnu.so $MYDIR
cp /g100/home/userexternal/aforment/.local/lib/python3.8/site-packages/sdf-1.0-py3.8.egg-info $MYDIR
```

check that you have installed sdf and sdf_helper by verifying that these commands don't give errors 
```
ipython3
import sdf
import sdf_helper
``` 