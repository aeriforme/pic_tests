---
title: Spack on laptop
updated: 2022-09-27 14:41:06Z
created: 2022-09-15 15:33:00Z
latitude: 45.46542190
longitude: 9.18592430
altitude: 0.0000
---

# What
software to install other software

# Documentation
https://spack.io/
https://spack.readthedocs.io/en/latest/
https://spack-tutorial.readthedocs.io/en/latest/

# Why
* useful to install dependencies 
* used on HPC systems 
* same package can be installed with different dependencies 
* environments 

# Install 
clone Spack repository 
```
git clone -c feature.manyFiles=true https://github.com/spack/spack.git spack
```
move to spack directory 
```
cd spack
```
load the setup script to add Spack to your path
```
source share/spack/setup-env.sh
```

# Usage
check your compilers 
```
spack compilers
```
you should see something like this 
```
==> Available compilers
-- gcc debian10-x86_64 ------------------------------------------
gcc@8.3.0  gcc@7.4.0  gcc@7.3.0

-- gcc debian11-x86_64 ------------------------------------------
gcc@10.2.1
```
choose your favorite compiler (maybe latest version) and always use that 

install the packages you need, for example this lets you install cmake version 3.23.3 using gcc version 10.2.1 
```
spack install cmake@3.23.3 %gcc@10.2.1
```

if you want to check out the available versions do this:  
```
spack versions cmake
``` 

another example
```
spack install hdf5 %gcc@10.2.1
```

hdf5 depends on an MPI implementation that, if not specified, Spack will choose for you 
you can specify which MPI implementation to use by explicitly adding the dependency with ^
this will install hdf5 using openmpi 
```
spack install hdf5 ^openmpi  %gcc@10.2.1
```
this will install hdf5 using mpich (openmpi and mpich are two different packages that provide an mpi implementation)
```
spack install hdf5 ^mpich  %gcc@10.2.1
```

you can enable a building option with + or disable it with ~
this explicitely says to install hdf5 enabling mpi using openmpi 
```
spack install hdf5+mpi ^openmpi  %gcc@10.2.1
```

if you want to disable mpi, i.e. install hdf5 without mpi support 
```
spack install hdf5~mpi
```

**Summary of the sigils** 
`%` = specify compiler 
`@` = specify version 
`^` = request dependency explicitely 
`+` = enable variant
`-` or `~` = disable variant 

# Do this
```
spack install cmake@3.23.3 %gcc@10.2.1
spack install openmpi %gcc@10.2.1
spack install hdf5+mpi ^openmpi %gcc@10.2.1
spack install adios2+mpi ^openmpi %gcc@10.2.1
spack install py-ipython
spack install openpmd-api +python -adios1 +adios2 -hdf5 +mpi
```

# Environments
create one environment for every pic code 
```
spack env create myspackwarpx
spack env create myspacksmilei
spack env create myspackepoch
```

activate the environment of the code you want to use, for example 
```
spack env activate -p myspackwarpx
```

to deactivate
```
despacktivate
```

in every environment install the required dependencies according to this list:
* WarpX
```
spack install cmake
spack install openmpi 
spack install adios2
spack install openpmd-api +python -adios1 +adios2 -hdf5 +mpi
```
* Smilei
```
spack install openmpi
spack install hdf5
spack install python
``` 
* EPOCH
```
spack install openmpi
``` 