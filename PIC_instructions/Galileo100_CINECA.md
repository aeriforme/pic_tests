# python
setup your own python environment
check out the instructions related to cineca

# WarpX
## Build
create a `warpx.profile` file in your `$HOME` with the following lines
```bash
module purge
module load gcc/10.2.0
module load cmake/3.21.4
module load profile/advanced
module load openmpi/4.1.1--gcc--10.2.0-cuda--11.1.0
module load zlib/1.2.11--gcc--10.2.0
module load adios/1.13.1--openmpi--4.1.1--gcc--10.2.0-cuda--11.1.0
module load boost/1.74.0--openmpi--4.1.1--gcc--10.2.0-cuda--11.1.0
```

source the file with 
```bash
source $HOME/warpx.profile
```

clone repository and move there 
```bash
git clone https://github.com/ECP-WarpX/WarpX.git warpx
cd warpx
```

then configure build
```bash
cmake -S . -B build
```

setup the compiling options (e.g. dimensionality, QED, table generation)
since not all the nodes in Galileo100 have GPUs, set `WarpX_COMPUTE = OMP` 
```bash
ccmake build
``` 

finalize building
```bash
cmake --build build -j 4
```

if successful the executable will be in `warpx/build/bin`

if you have problems with the pre-installed adios, you may want to try and install it by hand by following the instructions here: https://adios2.readthedocs.io/en/latest/setting_up/setting_up.html#install-from-source

```bash
git clone https://github.com/ornladios/ADIOS2.git ADIOS2
mkdir ADIOS2-build
cd ADIOS2-build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/ADIOS2-install ../ADIOS2
make -j 32
make install
```

## Run
move to `$CINECA_SCRACTH` and create a directory `MYDIR` there
then create a file `job.sh` with the following lines 
```bash
#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=31200MB 
#SBATCH -p g100_all_serial
#SBATCH --job-name=job_warpx
#SBATCH --err=warpx.err
#SBATCH --out=warpx.out
#SBATCH --account=pMI22_EneDa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mymail@mymail.com

export MYDIR=path/to/your/simulation/directory
cd $CINECA_SCRATCH/$MYDIR
source $HOME/warpx.profile

export OMP_SCHEDULE=dynamic
export OMP_NUM_THREADS=1

srun --cpu-bind=cores $HOME/warpx/build/bin/warpx.2d.MPI.OMP.DP.PDP.OPMD.QED.GENQEDTABLES input.txt > output.txt
```

submit your job with 
```bash
sbatch job.sh
```

## Post-process
you should already have installed the openpmd-api 
otherwise you can load your python environment and use pip

```bash
source $HOME/myenv/bin/activate
pip install openpmd-api
```

# Smilei
## Build
download the source code in your `$HOME`
```bash
git clone https://github.com/SmileiPIC/Smilei.git smilei
```

create a `smilei.profile` file in your `$HOME` with the following lines
```bash
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
```bash
source $HOME/smilei.profile
``` 

move to the right directory
```bash
cd smilei
``` 

compile 
```bash
make -j 8
```
if successful, you'll find the executable `smilei` in the current directory 
another executable `smilei_test` is generated to be run in test mode (e.g. check that the input file is ok)

## Run

prepare your job script `job.sh` and then submit it with `sbatch job.sh`

here is an example of what `job.sh` could contain if you want to run on the serial partition on 4 cores with 4 MPI tasks and 1 OMP thread per task 
```bash
#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=31200MB 
#SBATCH -p g100_all_serial
#SBATCH --job-name=job_smilei
#SBATCH --err=smilei.err
#SBATCH --out=smilei.out
#SBATCH --account=pMI21_EneDa_1
#SBATCH --qos=noQOS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mymail@mymail.com

cd $CINECA_SCRATCH/$MYDIR 

source $HOME/smilei.profile

export OMP_SCHEDULE=dynamic
export OMP_NUM_THREADS=1

srun --cpu-bind=cores -m block:block $HOME/smilei/smilei input.py > output.txt
```

here is an example of what `job.sh` could contain if you want to run on the serial partition on 4 cores with 1 MPI task and 1 OMP thread per core  
```bash
#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=31200MB 
#SBATCH -p g100_all_serial
#SBATCH --job-name=job_smilei
#SBATCH --err=smilei.err
#SBATCH --out=smilei.out
#SBATCH --account=pMI21_EneDa_1
#SBATCH --qos=noQOS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mymail@mymail.com

cd $CINECA_SCRATCH/$MYDIR 

source $HOME/smilei.profile

export OMP_SCHEDULE=dynamic
export OMP_NUM_THREADS=4

srun --cpu-bind=cores -m block:block $HOME/smilei/smilei input.py > output.txt
```


here it is assumed that: 
* you have created a directory `$MYDIR` in your scratch 
* in this directory you have created a valid `input.py` 
* you have the `smilei` directory and the `smilei.profile` file in your `$HOME`
* you have successfully compiled smilei in `$HOME/smilei/`

here are some tips for the parallelization from the smilei documentation: 
https://smileipic.github.io/Smilei/Understand/parallelization.html#practical-setup

:warning: warning: performances are extremely sensitive to the number of patches you use (which has to be a power of 2 in each direction) &rarr; you may have to redefine the number of cells in your simulation to optimize the number of patches you can use (often, the more the better)

## Post-process
load `smilei.profile` and your python environment, then move to the source directory
```bash
source $HOME/smilei.profile
source $HOME/myenv/bin/activate
cd $HOME/smilei
```

build the post-processing python package `happi`
```bash
make happi
```
a message like this should appear
```bash
Installing /g100/home/userexternal/<username>/.local/lib/python3.8/site-packages/smilei.pth
```

copy the smilei.pth path in the correct directory of your python environment 
```bash
cp /g100/home/userexternal/<username>/.local/lib/python3.8/site-packages/smilei.pth $HOME/myenv/lib/python3.8/site-packages
```

now you should be able to `import happi` in a python shell or script &rarr; you can use your plot files to post-process your smilei data as you would do on your computer (provided that you source your python environment first)

# EPOCH
## Build
clone the repository
```bash
git clone --recursive https://github.com/Warwick-Plasma/epoch.git
```

create a `epoch.profile` file in your `$HOME` with the following lines
```bash
module purge
module load intel/oneapi-2022--binary
module load intelmpi/oneapi-2022--binary
module load python/3.8.12--intel--2021.4.0
```

then move to the target directory, e.g.
```bash
cd epoch/epoch2d
```

and compile by typing
```bash
source $HOME/epoch.profile
make -j8 COMPILER=intel
```

the executable should be generated in `epoch/epoch2d/bin`

## Run
prepare your job script `job.sh` and then submit it with `sbatch job.sh`

here is an example of what `job.sh` could contain 
```bash
#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=31200MB 
#SBATCH -p g100_all_serial
#SBATCH --job-name=job_epoch
#SBATCH --err=epoch.err
#SBATCH --out=epoch.out
#SBATCH --account=pMI21_EneDa_1
#SBATCH --qos=noQOS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mymail@mymail.com

cd $CINECA_SCRATCH/$MYDIR
source $HOME/epoch.profile
echo Data | srun $HOME/epoch/epoch2d/bin/epoch2d > output.txt
```

here it is assumed that: 
* you have created a directory `$MYDIR` in your scratch 
* in this directory you have created a `Data` sub-directory 
* you have a valid `input.deck` inside the `Data` directory 
* you have the `epoch` directory and the `epoch.profile` file in your `$HOME`
* you have successfully compiled epoch2d

## Post-process 
source your python3 environment 
move to any epoch source directory (e.g. `epoch/epoch2d`) and type
```bash
make sdfutils
```
a message similar to this should appear
```bash
writing list of installed files to 'pybuild/files.txt'
```
you have to copy all the files listed in `/epoch/SDF/utilities/pybuild/files.txt` in the correct directory of your python environment `$HOME/myenv/lib/python3.8/site-packages/`
this should work
```bash
export MYDIR=$HOME/myenv/lib/python3.8/site-packages/
cp /g100/home/userexternal/aforment/.local/lib/python3.8/site-packages/sdf_legacy.py $MYDIR
cp -r /g100/home/userexternal/aforment/.local/lib/python3.8/site-packages/sdf_helper $MYDIR/sdf_helper
cp /g100/home/userexternal/aforment/.local/lib/python3.8/site-packages/__pycache__/sdf_legacy.cpython-38.pyc $MYDIR/__pycache__
cp /g100/home/userexternal/aforment/.local/lib/python3.8/site-packages/sdf.cpython-38-x86_64-linux-gnu.so $MYDIR
cp /g100/home/userexternal/aforment/.local/lib/python3.8/site-packages/sdf-1.0-py3.8.egg-info $MYDIR
```

check that you have installed sdf and sdf_helper by verifying that these commands don't give errors 
```bash
source $HOME/myenv/bin/activate
ipython3
import sdf
import sdf_helper
``` 

now you can use the python scripts to post-process your epoch data as you would do on your computer (provided that you source your python environment)