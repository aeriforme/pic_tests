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
clone repository and move there 
```
git clone https://github.com/ECP-WarpX/WarpX.git warpx
cd warpx
```

## Run

## Post-process


# Smilei
## Build
download the source code in your `$HOME`
```
git clone https://github.com/SmileiPIC/Smilei.git smilei
```

create a `smilei.profile` file in your `$HOME` with the following lines
```
module purge 
module load intel/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary
module load libszip/2.1.1--gcc--10.2.0
module load zlib/1.2.11--gcc--10.2.0
module load hdf5/1.10.7--intelmpi--oneapi-2021--binary
module load anaconda3/2020.07--gcc--8.3.1
export PYTHONEXE=python3
export HDF5_ROOT_DIR=/cineca/prod/opt/libraries/hdf5/1.10.7/intelmpi--oneapi-2021--binary/
export SMILEICXX=mpiicpc
```

and source it
```
source $HOME/smilei.profile
``` 

move to the right directory
```
cd smilei
``` 

compile 
```
make -j 8
```
if successful, you'll find the executable `smilei` in the current directory 
another executable `smilei_test` is generated to be run in test mode (e.g. check that the input file is ok)

## Run

## Post-process
load `smilei.profile` and your python environment, then move to the source directory
```
source $HOME/smilei.profile
source $HOME/myenv/bin/activate
cd $HOME/smilei
```

build the post-processing python package `happi`
```
make happi
```
a message like this should appear
```
Installing /g100/home/userexternal/<username>/.local/lib/python3.8/site-packages/smilei.pth
```

copy the smilei.pth path in the correct directory of your python environment 
```
cp /g100/home/userexternal/<username>/.local/lib/python3.8/site-packages/smilei.pth $HOME/myenv/lib/python3.8/site-packages
```

now you should be able to `import happi` in a python shell or script 

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

## Run
create a directory in $CINECA_SCRATCH
cd $CINECA_SCRATCH
mkdir epoch_tests

copy an input file in this directory 

## Post-process 
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