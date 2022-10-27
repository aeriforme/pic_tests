# Documentation 
* general: https://wiki.u-gov.it/confluence/display/SCAIUS/HPC+at+CINECA%3A+User+Documentation 
* Galileo100: https://wiki.u-gov.it/confluence/display/SCAIUS/UG3.3%3A+GALILEO100+UserGuide
* Marconi100: https://wiki.u-gov.it/confluence/display/SCAIUS/UG3.2%3A+MARCONI100+UserGuide
* SLURM scheduler: https://slurm.schedmd.com/documentation.html

(and soon the new supercomputer Leonardo...)

# Software
you can load (e.g. set your environment variables appropriately) via the module command 
```text
Command                Action
---------------------------------------------------------------------------------------------- 
module avail .................. show the available modules on the machine 
module load <appl> ............ load the module <appl> in the current shell session, preparing the environment for the application. 
module load autoload <appl> ... load the module <appl> and all dependencies in the current shell session
module help <appl> ............ show specific information and basic help on the application 
module list ................... show the modules currently loaded on the shell session 
module purge .................. unload all the loaded modules 
module unload <appl>........... unload a  specific module
modmap -m <appl>................find available <appl> 
```

if you need some software that is not pre-installed on the cineca machines, you should try installing it with spack (which is installed at cineca and you can use it after loading it: `module load spack`)

if spack doesn't provide that software or if it doesn't work, you will likely have to install it from source (could be painful)

# File transfer

to transfer files from a supercomputer to your local computer or viceversa you can use `scp` or `rsync`

to copy from remote to local, from your laptop do 
```bash
scp <username>@login.m100.cineca.it:<path_2_file_u_want_2_download> <path_to_local_directory>
```
or
```bash
rsync --update --progress <username>@login.m100.cineca.it:<path_2_file_u_want_2_download>  <path_to_local_directory>
```

to copy from local to remote, from your laptopo do
```bash
scp <path_2_file_u_want_2_upload> <username>@login.m100.cineca.it:<your_target_directory> 
```
or
```bash
rsync --update --progress <path_2_file_u_want_2_upload> <username>@login.m100.cineca.it:<your_target_directory>
```

# ssh
instructions to connect to the supercomputers via ssh can be found here:
https://www.hpc.cineca.it/content/how-connect-public-key

here a brief recap:
type in  your terminal the following command to generate a ssh key
```bash
ssh-keygen
```
press enter when prompted 
then you will find two files in your `$HOME/.ssh` directory:`id_rsa` and `id_rsa.pub` (check this by doing `ls ~/.ssh`) 

now copy `id_rsa.pub` on the remote server, e.g. on Galileo100
```bash
scp ~/.ssh/id_rsa.pub <username>@login.g100.cineca.it:/g100/home/userexternal/<username>
```

and add it to the authorized keys with 
```bash
cd $HOME
mkdir -p .ssh 
cat id_dsa.pub >> $HOME/.ssh/authorized_keys
```

if successful, now you're not asked to type the password when you connect

you can also create aliases in your `~/.bashrc` to connect more easily, for example: 
```bash
alias galileo100='ssh -X <username>@login.g100.cineca.it'
alias marconi100='ssh -X <username>@login.m100.cineca.it'
```
then you can connect to Galileo100 by just typing `galileo100` on your terminal 

# Filesystems
more info here: https://wiki.u-gov.it/confluence/display/SCAIUS/UG2.5%3A+Data+storage+and+FileSystems

### `$HOME`
* it's where you land when you log in 
* permanent: files won't be removed 
* use this space to store your sources and binaries and other lightweight material
* compile your codes here, but do not run simulations here (unless extremely tiny, e.g. just to check that a binary works with a given input)
* 50 GB quota
* it's backed up 
* user-specific (and private)

### `$CINECA_SCRATCH`
* move there by `cd $CINECA_SCRATCH`
* user-specific 
* private, but you can change the permissions of this directory (at your own risk) if you want to easily exchange files with some other user by being able to read and browse their scratch. by default  the permissions of your directory in the scratch should be `drwx------` (this is the output of  `ls -ld`). you can change this to `drwxr-xr-x`  with `chmod`:
`chmod 755 <your_directory>`: now everybody should be able to read and browse `<your_directory>`

*  temporary: there is a periodic (daily) cleaning procedure: files that have remained untouched for more than 40 days are removed and a file in your home is generated listing all the deleted files on a given day (if any). to avoid losing your files unintentionally, you can touch all files and subdirectories with this command `find . -exec touch {} \;` (this however is frowned upon!)
*  "no disk quota" but a temporary quota of 20 TB could be imposed if the cluster reaches 88% of occupancy 

###`$WORK`
* project-specific (you share this with your colleagues, so talk to them to decide how to use it) 
* collaborative space 
* quota = 1 TB


to get some stats of your disk storage:
```bash
cindata
```

to get stats on your budget situation 
```bash
saldo -b
saldo -r
``` 


# Jobs
## batch 
to run your simulations you have to submit jobs to the scheduler
the scheduler will put your job in queue, assign the computational resources required, launch it, etc.

```bash
#!/bin/bash 

# slurm directives
#SBATCH --time=04:00:00 # walltime 
#SBATCH --nodes=1 # number of nodes 
#SBATCH --ntasks-per-node=4 # number of cores
#SBATCH --cpus-per-task=1 # number of mpi tasks/core
#SBATCH --mem=31200MB # memory/core
#SBATCH -p g100_all_serial # partition
#SBATCH --qos=noQOS # quality of service
#SBATCH --job-name=myname 
#SBATCH --err=myjob.err # standard error
#SBATCH --out=myjob.out # standard output
#SBATCH --account=pMI21_EneDa_1 # account / budget
#SBATCH --mail-type=ALL # type of alerts
#SBATCH --mail-user=mymail@mymail.com # email adress

# move to the directory where you want to run the code
cd /g100_scratch/userexternal/aforment/tests/pic_tests/2d_laser_plasma/epoch

# source the required software, e.g. for epoch 
source epoch.profile

# run the code e.g. with srun 
echo Data | srun ./epoch2d > output.txt
```

to submit your job 
```bash
sbatch job.sh
``` 
to check the status of your job 
```bash
squeue -u <username>
```
you will see some stats of your job: the ID assigned by the scheduler, the status, the running time, etc.  

to get some other stats of your job do this
```bash
scontrol show job <job_id>
```

if you want to cancel your
```bash
scancel <jobid>
```

## interactive
 
you can also ask for computational resources interactively: this means that you will be able to type in commands on a compute node 

for example:
```bash
srun -n1 --time=04:00:00 -A pMI21_EneDa_1 -p g100_all_serial --ntasks-per-node=1 --cpus-per-task=1 --pty /bin/bash
```

when you close your terminal session, the allocation will end 

# python environment

create your own python environment in your $HOME
```bash
module load python/3.8.12--gcc--10.2.0
python3 -m venv myenv
source $HOME/myenv/bin/activate
```

check that you can use your python environment with
```bash
which python3
which pip3
```

now install the required software
```bash
pip install numpy 
pip install scipy 
pip install ipython 
pip install matplotlib 
pip install openpmd-api
``` 
