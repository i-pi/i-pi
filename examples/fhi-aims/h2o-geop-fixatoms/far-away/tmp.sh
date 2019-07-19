#!/bin/bash

for mode in bfgs lbfgs bfgstrm cg sd
do
    for i in `seq 0 1`
    do
        mv aims.${i}-fixed.${mode}.faraway.out  aims.${i}-fixed.${mode}.out
        mv sim-0-fixed.${mode}.faraway.pos_0.xyz sim.${i}-fixed.${mode}.pos_0.xyz
    done
done
