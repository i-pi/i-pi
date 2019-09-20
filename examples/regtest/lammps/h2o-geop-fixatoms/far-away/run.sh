#!/bin/bash

. ../../../../../env.sh
IPI_EXE='python2 -u ../../../../../bin/i-pi'
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

# 0 fixed atoms
sed -e "s/PREFIX/sim-0-fixed.${mode}/g" input.template.xml > input.xml
sed -i -e "s/^.*FIXATOMS.*$//g" input.xml
sed -i -e "s/MODE/${mode}/g" input.xml
sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
subdir="${mode}-0-fixed"
mkdir $subdir
cp {input.xml,init.far-away.xyz,data.water,in.lmp} $subdir/

# 1 fixed atom
sed -e "s/PREFIX/sim-1-fixed.${mode}/g" input.template.xml > input.xml
sed -i -e "s/FIXATOMS/3/g" input.xml
sed -i -e "s/MODE/${mode}/g" input.xml
sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
subdir="${mode}-1-fixed"
mkdir $subdir
cp {input.xml,init.far-away.xyz,data.water,in.lmp} $subdir/

rm input.xml
