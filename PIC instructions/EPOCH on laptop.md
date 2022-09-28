# Documentation
https://epochpic.github.io/

# Build 
clone the repository 
```
git clone --recursive https://github.com/Warwick-Plasma/epoch.git
```

move to the directory 
```
cd epoch
```

and move again to choose in which dimensionality you want to compile the code (in EPOCH there's different source codes for 1D, 2D and 3D geometry)
for example to compile the 2D version of EPOCH do this 
```
cd epoch2d
```

source Spack if you don't do it in your `bash.rc` 
```
source $SPACK_DIR/share/spack/setup-env.sh
```

load the epoch enviroment 
```
spack env activate -p myspackepoch
```

compile
```
make COMPILER=gfortran
```
once this command is done you will have your epoch2d executable in `epoch/epoch2d/bin`

you can activate a number of compile-time options by modifying the `Makefile` 
you can clean the building process with `make cleanall` (do it everytime you re-compile, just to be sure)

and build the python library to read epoch output (called sdf)
```
make sdfutils
```

in a python script or shell you should now be able to use 
```
import sdf 
import sdf_helper as sh
``` 


# Run 
make a directory for your epoch tests `tests_epoch`
```
mkdir tests_epoch
``` 

make a subdirectory `Data`
```
cd tests_epoch 
mkdir Data
``` 

copy in `tests_epoch` the epoch executable 
copy in  `tests_epoch/Data`  an input file `input.deck`

move to `tests_epoch` and run on 4 cores, e.g. 
```
echo Data | mpirun -np 4 ./epoch2d
```
