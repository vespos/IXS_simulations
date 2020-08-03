#!/bin/bash

# seq FIRST STEP LAST
corrs=($(seq -0.005 0.001 0.015))
# corrs=($(seq 0 0.001 0.001))

FILE=./dat/optim_test.dat
if [ -f "$FILE" ]; then
    rm $FILE
fi

for corr in "${corrs[@]}"
do
    echo Simulating pitch correction: $corr
    python ./spectro_9keV_pitch_optim.py -a $corr
done