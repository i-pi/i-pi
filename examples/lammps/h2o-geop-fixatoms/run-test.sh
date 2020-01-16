#!/bin/bash

. ~/soft/i-pi-mahrossi/env.sh
IPI_EXE='python3 -u /Users/rossi/Codes/i-pi/bin/i-pi'
#IPI_EXE='python2 -u /home/fidanyan/soft/i-pi-cosmo/bin/i-pi'
LMP_EXE=lmp_serial

mode=$1
if [ $mode -eq '']
then
    echo -e "ERROR: 1 argument needed:\n\tmode - bfgs, bfgstrm, lbfgs, cg or sd."
    exit
fi

rm EXIT

# 2 fixed atoms
sed -e "s/PREFIX/sim-2-fixed.${mode}/g" input.template.xml > input.xml
sed -i -e "s/FIXATOMS/1,2/g" input.xml
sed -i -e "s/MODE/${mode}/g" input.xml
$IPI_EXE input.xml |tee log.2-fixed.${mode}.ipi &

echo
echo "2..."
sleep 1
echo "1..."
sleep 1
$LMP_EXE -in in.lmp -log log.2-fixed.${mode}.lmp &
wait

# For a comparison, here you can run unconstrained optimization
# 0 fixed atoms
sed -e "s/PREFIX/sim-0-fixed.${mode}/g" input.template.xml > input.xml
sed -i -e "/FIXATOMS/d" input.xml
sed -i -e "s/MODE/${mode}/g" input.xml
$IPI_EXE input.xml |tee log.0-fixed.${mode}.ipi &

echo
echo "2..."
sleep 1; echo "1..."; sleep 1
$LMP_EXE -in in.lmp -log log.0-fixed.${mode}.lmp &
wait
