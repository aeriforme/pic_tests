---
title: EPOCH on laptop
updated: 2022-09-27 14:26:22Z
created: 2022-09-15 16:36:07Z
latitude: 45.46542190
longitude: 9.18592430
altitude: 0.0000
---

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

and build the python library to read EPOCH output (called sdf)
```
make sdfutils
```


# Run 
run on 4 cores 
```
echo Data | mpirun -np 4 ./epoch2d
```

# Analyze
in a python script you can use 
```
import sdf 
import sdf_helper as sh
``` 
