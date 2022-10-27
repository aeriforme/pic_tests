You can choose to install the main dependencies 
* from source
* via the package manager of your system
* using Spack (another package management tool) 

# Via apt or source 
on a Debian or Ubuntu machine install the following software 
```
sudo apt-get install git python3-h5py python3-ipython python3-pint python3-sphinx python3-matplotlib python3-dev python3-numpy python3-pip build-essential gcc libhdf5-openmpi-dev
```

you may need to install openmpi and/or hdf5 from source 
to do that, folllow these instructions: https://smileipic.github.io/Smilei/Use/install_linux.html#troubleshooting

# Spack
### What
software to install other software that is useful because 
*  works on different platforms and environments like Linux machines, macOS and HPC systems (and Windows is work in progress)
* the same package can be installed using different dependencies, compilers and options without conflicts 
* environments can make the workflow conveniently modular 

### Documentation
https://spack.io/
https://spack.readthedocs.io/en/latest/
https://spack-tutorial.readthedocs.io/en/latest/

### Install 
clone Spack repository 
```
git clone -c feature.manyFiles=true https://github.com/spack/spack.git spack
```
move to Spack directory 
```
cd spack
```
load the setup script to add Spack to your path
```
source share/spack/setup-env.sh
```

you can add the source line in your `.bashrc` to avoid having to do it every time 

### Usage
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

### What to install
Install the following packages 
```
spack install cmake@3.23.3 %gcc@10.2.1
spack install openmpi %gcc@10.2.1
spack install hdf5+mpi ^openmpi %gcc@10.2.1
spack install adios2+mpi ^openmpi %gcc@10.2.1
spack install openpmd-api+mpi+python+adios2-adios1 %gcc@10.2.1  
spack install py-pip
```

if you want to use 
pip install matplotlib numpy scipy ipython pint sphinx h5py openpmd-api 


cmake and openmpi (plus git and a C++ compiler) are dependencies common to the 3 PIC codes (WarpX, Smilei and EPOCH)
hdf5 is mandatory for Smilei and optional for WarpX
adios2 is optional for WarpX

note: WarpX itself also exists as a Spack package 
```
spack info warpx
``` 
to install, do for example
```
spack install warpx +sensei dims=2 compute=omp
``` 

### Environments
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