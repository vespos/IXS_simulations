#!/bin/bash

# seq FIRST STEP LAST
pitchs=($(seq 0.5 0.1 1.5))

FILE=./dat/ML.dat
if [ -f "$FILE" ]; then
    rm $FILE
fi

for pitch in "${pitchs[@]}"
do
    echo Simulating pitch: $pitch
    python ./multilayer_reflectivity.py -a $pitch -f $FILE
done