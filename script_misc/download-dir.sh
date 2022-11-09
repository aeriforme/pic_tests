#!/bin/bash

name=BREMBELLA
source=aforment@login.g100.cineca.it:/g100_scratch/userexternal/aforment/${name}
echo ${source}

dest=${name}
while [ 1 ]
do
    rsync -vv --partial --progress --inplace --append --compress --recursive $source $dest
    if [ "$?" = "0" ] ; then
        echo "rsync completed normally"
        exit
    else
        echo "Rsync failure. Backing off and retrying..."
        sleep 180
    fi
done

