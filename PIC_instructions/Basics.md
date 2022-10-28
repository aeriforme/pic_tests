# Goal of these instructions
help one build and run 3 particle-in-cell codes (WarpX, Smilei and EPOCH) on a laptop with linux or on Cineca supercomputers

refer to this github repository https://github.com/aeriforme/pic_tests for examples and material

# Laptop
## General dependencies
fist, install on you laptop the following fundamental software
```bash
sudo apt-get install git python3-h5py python3-ipython python3-pint python3-sphinx python3-matplotlib python3-dev python3-numpy python3-scipy python3-pip build-essential gcc libhdf5-openmpi-dev 
```
then, follow the related instructions for each case

to run efficiently in parallel on your machine you need to know your architecture
for example, you can find out the number of threads per core and cores per socket of your machine with the command `lscpu`

# Cineca
you must have a cineca account

to access Galileo100:
```bash
ssh <username>@login.g100.cineca.it
```
similar for Marconi100
```bash
ssh <username>@login.m100.cineca.it
```
then follow the related instructions 

# Git 
there's plenty of tutorials online, please consider this as just a summary/reminder of the main git commands 

if you want to contribute to the `pic_tests` repo, first fork it in your github account, then clone it 
```bash
git clone https://github.com/<your_git_username>/pic_tests
```

to incorporate changes from a remote repository into the current branch 
```bash
git pull
```

to prepare the content staged for the next commit
```bash
git add <files_to_add>
```

create a new commit, i.e. record the latest changes of the source code to your repository
```bash
git commit -a -m "useful message"
```

actually update remote repository
```bash
git push
```

if you want to merge your changes in the main repository, go on github and make a pull request

if the pull request is accepted and the code is merged, then go back to your repo and synch your fork + do git pull from terminal to update your code 
