# Documentation
https://smileipic.github.io/Smilei

# Build 
download the source code
```
git clone https://github.com/SmileiPIC/Smilei.git smilei
```

move to the right directory
```
cd smilei
``` 

source spack if you don't do it in your `bash.rc` 
```
source $SPACK_DIR/share/spack/setup-env.sh
```

load the smilei enviroment 
```
spack env activate -p myspacksmilei
```

set the environment variable `HDF5_ROOT_DIR` to the path where the correct installation of hdf5 is  
you can find the path in the output of the following command 
```
spack find --path hdf5+mpi ^openmpi %gcc@10.2.1
```
for example in my case
```
export HDF5_ROOT_DIR =$HOME/src/spack/opt/spack/linux-debian11-skylake/gcc-10.2.1/hdf5-1.12.2-aqc4gaolbexiqq7wvne2bhewi6fjxet3
export LD_LIBRARY_PATH=${HDF5_ROOT_DIR}/lib:${LD_LIBRARY_PATH}
```

set also this other environment variable 
```
export PYTHONEXE=python3
``` 

you can set the environment variables in your `.bashrc` file

then compile
```
make -j 2
``` 

build the python module to manage output 
```
make happi
```

in a python script you can now do
```
import happi
```

# Run 
copy the executable (smilei) and an input file (input.py) in a directory (mydir) and move there 

set the number of threads per core depending on the machine 
for example, if Thread(s) per core = 2 (in the output of the `lscpu` command), then
```
export OMP_NUM_THREADS=2
```

and if you have Core(s) per socket = 4, you can run on the 4 cores like this 
```
mpirun -np 4 ./smilei input.py
```