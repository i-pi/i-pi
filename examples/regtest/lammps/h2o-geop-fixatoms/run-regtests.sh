#!/bin/bash

. ../../../../env.sh
IPI_EXE='python2 -u ../../../../bin/i-pi'
LMP_EXE=lmp_mpi

mode=$1
if [ $mode -eq '']
then
    echo -e "ERROR: 1 argument needed:\n\tmode - bfgs, bfgstrm, lbfgs, cg or sd."
    exit
fi

rm EXIT

gitbranch=`git branch|grep '*'|awk '{print $2}'`
echo "Current i-PI git branch is $gitbranch."

#mpirun -np 4 $LMP_EXE -in in.lmp -log log.0-fixed.${mode}.lmp &

# 0 fixed atoms
sed -e "s/PREFIX/sim-0-fixed.${mode}/g" input.template.xml > input.xml
sed -i -e "s/^.*FIXATOMS.*$//g" input.xml
sed -i -e "s/MODE/${mode}/g" input.xml
sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
subdir="${mode}-0-fixed"
mkdir $subdir
cp {input.xml,init.xyz,data.water,in.lmp} $subdir/

# 1 fixed atom
sed -e "s/PREFIX/sim-1-fixed.${mode}/g" input.template.xml > input.xml
sed -i -e "s/FIXATOMS/1/g" input.xml
sed -i -e "s/MODE/${mode}/g" input.xml
sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
subdir="${mode}-1-fixed"
mkdir $subdir
cp {input.xml,init.xyz,data.water,in.lmp} $subdir/

# 2 fixed atoms
sed -e "s/PREFIX/sim-2-fixed.${mode}/g" input.template.xml > input.xml
sed -i -e "s/FIXATOMS/1,2/g" input.xml
sed -i -e "s/MODE/${mode}/g" input.xml
sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
subdir="${mode}-2-fixed"
mkdir $subdir
cp {input.xml,init.xyz,data.water,in.lmp} $subdir/


# 3 fixed atoms (all)
sed -e "s/PREFIX/sim-3-fixed.${mode}/g" input.template.xml > input.xml
sed -i -e "s/FIXATOMS/0,1,2/g" input.xml
sed -i -e "s/MODE/${mode}/g" input.xml
sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
subdir="${mode}-3-fixed"
mkdir $subdir
cp {input.xml,init.xyz,data.water,in.lmp} $subdir/

cd far-away
    ./run.sh ${mode}
cd ..

rm input.xml
../../../../tools/py/regtest.py --create-reference
