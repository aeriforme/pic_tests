# Goal of these instructions
help one building and running 3 particle-in-cell codes (WarpX, Smilei and EPOCH) on a laptop with linux or on Cineca supercomputers

refer to this github repository https://github.com/aeriforme/pic_tests for examples and material

# Laptop
## General dependencies
fist, install on you laptop the following fundamental software
```
sudo apt-get install git python3-h5py python3-ipython python3-pint python3-sphinx python3-matplotlib python3-dev python3-numpy python3-pip build-essential gcc libhdf5-openmpi-dev
```
then, follow the related instructions

# Cineca
you must have a cineca account. 
to access Galileo100:
```
ssh <username>@login.g100.cineca.it
```
similar for Marconi100
```
ssh <username>@login.m100.cineca.it
```
then follow the related instructions 

to transfer files from a supercomputer to your local computer or viceversa you can use `scp` or `rsync`

to copy from remote to local, from your laptop do 
```
scp <username>@login.m100.cineca.it:<your_directory>/<your_file> <path_to_local_directory>
```
or
```
rsync --update --progress <username>@login.m100.cineca.it:<your_directory>/<your_file> <path_to_local_directory>
```

to copy from local to remote, from your laptopo do
```
scp <path_to_local_file> <username>@login.m100.cineca.it:<your_target_directory> 
```
or
```
rsync --update --progress <path_to_local_file> <username>@login.m100.cineca.it:<your_target_directory>
```
