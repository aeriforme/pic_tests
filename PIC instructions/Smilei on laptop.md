---
title: Smilei on laptop
updated: 2022-09-27 14:21:03Z
created: 2022-09-15 16:58:54Z
latitude: 45.46542190
longitude: 9.18592430
altitude: 0.0000
---

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

source Spack if you don't do it in your `bash.rc` 
```
source $SPACK_DIR/share/spack/setup-env.sh
```

load the Smilei enviroment 
```
spack env activate -p myspacksmilei
```

compile
```
make -j 2
``` 

build the python module to manage output 
```
make happi
```

# Run 
copy the executable (smilei) and an input file (input.py) in a directory (mydir) and move there 
set the number of threads per core depending on the machine 
```
export OMP_NUM_THREADS=2
```
run on 4 cores 
```
mpirun -np 4 ./smilei input.py
```

# Analyze
in a python script you can now do
```
import happi
```
