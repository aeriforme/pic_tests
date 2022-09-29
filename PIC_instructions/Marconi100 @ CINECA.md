# Build WarpX
clone repository and move there 
```
git clone https://github.com/ECP-WarpX/WarpX.git warpx
cd warpx
```

create a `warpx.profile` file to setup the environment by loading the modules that are dependencies  
put this file in your `$HOME` directory

but first, install adios2 using spack, which is provided by cineca
do this: 
```
module purge
module load gnu/8.4.0
module load cmake/3.20.0
module load spectrum_mpi/10.3.1--binary
module load cuda/11.3
module load python/3.8.2
module load boost/1.72.0--spectrum_mpi--10.3.1--binary
module load zlib/1.2.11--gnu--8.4.0
spack install adios2+mpi %gcc@8.4.0
```
copy the path of the directory where adios2 is installed (you can find it at the end of spack's output)

then copy everything into the `warpx.profile` file 
```
module purge

module load gnu/8.4.0
module load cmake/3.20.0
module load spectrum_mpi/10.3.1--binary
module load cuda/11.3
module load python/3.8.2
module load boost/1.72.0--spectrum_mpi--10.3.1--binary
module load zlib/1.2.11--gnu--8.4.0

#spack install adios2+mpi %gcc@8.4.0
export ADIOS2_DIR=/m100_work/pMI21_EneDa_0/aforment/spack-0.14/install/linux-rhel8-power9le/gcc-8.4.0/adios2-2.5.0-5buwthqfrizkisr4xsba6g6bmeb4hvau
```

then configure build
```
cmake -S . -B build
```

setup the compiling options (e.g. dimensionality, QED)
since Marconi100 has 4 GPUs on each node, set `WarpX_COMPUTE = CUDA` 
```
ccmake build
``` 

finalize building
```
cmake --build build -j 4
```

if successful the executable will be in `warpx/build/bin`

# Run WarpX
create a directory in `$CINECA_SCRATCH` where you put an input file and a warpx executable, here `warpx.2d.MPI.CUDA.DP.PDP.OPMD.QED.GENQEDTABLE`
create a job.x file with the following
```
#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4 #32 
#SBATCH --ntasks-per-socket=2
#SBATCH --gpus-per-node=4
#SBATCH --mem=61500 #7100 #246000
#SBATCH -p m100_usr_prod
#SBATCH --qos=m100_qos_dbg
#SBATCH --job-name=warpx
#SBATCH --err=job.err
#SBATCH --out=job.out
#SBATCH --account=pMI21_eneDa_0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arianna.formenti@polimi.it

export GPUS_PER_SOCKET=2
export GPUS_PER_NODE=4

# optimize CUDA compilation for NVIDIA V100 (7.0) or for A100 (8.0)
export AMREX_CUDA_ARCH=7.0

source ~$HOME/warpx.profile
cd $CINECA_SCRATCH/test

mpirun -gpu --rank-by core -n 4 $HOME/newwarpx/warpx/build/bin/warpx.2d.MPI.CUDA.DP.PDP.OPMD.QED.GENQEDTABLES input.txt warpx.numprocs=1 4 > output.txt
```

submit the job with 
```
sbatch job.x
```

check your jobs with 
```
squeue -u <user>
scontrol show job <job_id>
``` 

cancel your job with 
```
scancel <job_id>
```

# Post-process 
create your python3 environment in your `$HOME` directory 
```
module load python/3.8.2
python3 -m venv myenv
source myenv/bin/activate
```

modify the script `myenv/bin/activate`  by adding on top the following lines 
```
module load gnu/8.4.0
module load cmake/3.20.0
module load spectrum_mpi/10.3.1--binary
module load cuda/11.3
module load python/3.8.2
module load openblas/0.3.9--gnu--8.4.0
module load numpy/1.19.4--python--3.8.2
module load blas/3.8.0--gnu--8.4.0
module load lapack/3.9.0--gnu--8.4.0
module load profile/deeplrn
module load scipy/1.5.1--python--3.8.2
```

and then re-sources it 
```
source myenv/bin/activate
```

install ipython, just for convenience  
```
pip3 install ipython
``` 
and the openpmd api
```
pip3 install openpmd-api
``` 


## vtk 
if you need vtk, you have to install it from source 
follow the instructions here 
https://gitlab.kitware.com/vtk/vtk/-/blob/master/Documentation/dev/build.md

```
git clone --recursive https://gitlab.kitware.com/vtk/vtk.git
```

load the dependencies 
```
module purge
module load gnu/8.4.0
module load cmake/3.20.0
module load spectrum_mpi/10.3.1--binary
module load cuda/11.3
module load python/3.8.2
```

```
cd vtk 
mkdir -p build
cd build
ccmake .. # select e.g MPI, CUDA, PYTHON
cmake -DVTK_WHEEL_BUILD=ON -DVTK_WRAP_PYTHON=ON ../
cmake --build . -j8
```

```
source myenv/bin/activate
python3 setup.py bdist_whee # this generates a .whl file in /build/dist
pip3 install wheel build/dist/file.whl
```
