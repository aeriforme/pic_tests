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

you can add the source line in your `.bashrc` to avoid having to do it every time 

# Usage
add your compilers 
```
spack compiler find
``` 

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
choose your favorite compiler among those in your list (maybe latest version) and always use that 

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
spack install openpmd-api +python -adios1 +adios2 -hdf5 +mpi %gcc@10.2.1
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
or, its abbreviation, 
```
spacktivate -p myspackwarpx
```

to list all the environments you have in your system do this 
```
spack env list
``` 

to check what environment you are in, do 
```
spack env status
``` 

to deactivate
```
despacktivate
```

in every environment install the required dependencies according to this list:
* WarpX
```
spack install cmake@3.23.3 %gcc@10.2.1
spack install openmpi %gcc@10.2.1 
spack install adios2+mpi ^openmpi %gcc@10.2.1
spack install openpmd-api +python -adios1 +adios2 -hdf5 +mpi %gcc@10.2.1
```
* Smilei
```
spack install openmpi %gcc@10.2.1 
spack install hdf5+mpi ^openmpi %gcc@10.2.1
``` 
* EPOCH
```
spack install openmpi %gcc@10.2.1
``` 

note that doing `spack install` within an environment does not re-install the software from the beginning, it associates a specific install to your environment 