#!/bin/bash
# Helper script to create all test cases if I need to make some changes

for mode in bfgs bfgstrm lbfgs cg sd
do
    # 0 fixed atoms
    sed -e "s/PREFIX/sim-0fixed.${mode}/g" analysis-scripts/input.template.xml > input.xml
    sed -i -e "s/^.*FIXATOMS.*$//g" input.xml
    sed -i -e "s/MODE/${mode}/g" input.xml
    sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
    subdir="${mode}-0fixed"
    echo $subdir
    if [ ! -d $subdir ]; then mkdir $subdir; fi
#    cp {input.xml,init.xyz,data.water,in.lmp} $subdir/
    cp input.xml $subdir/

    # 1 fixed atom
    sed -e "s/PREFIX/sim-1fixed.${mode}/g" analysis-scripts/input.template.xml > input.xml
    sed -i -e "s/FIXATOMS/1/g" input.xml
    sed -i -e "s/MODE/${mode}/g" input.xml
    sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
    subdir="${mode}-1fixed"
    echo $subdir
    if [ ! -d $subdir ]; then mkdir $subdir; fi
#    cp {input.xml,init.xyz,data.water,in.lmp} $subdir/
    cp input.xml $subdir/

    # 2 fixed atoms
    sed -e "s/PREFIX/sim-2fixed.${mode}/g" analysis-scripts/input.template.xml > input.xml
    sed -i -e "s/FIXATOMS/1,2/g" input.xml
    sed -i -e "s/MODE/${mode}/g" input.xml
    sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
    subdir="${mode}-2fixed"
    echo $subdir
    if [ ! -d $subdir ]; then mkdir $subdir; fi
#    cp {input.xml,init.xyz,data.water,in.lmp} $subdir/
    cp input.xml $subdir/

    # Far away, 0 fixed atoms
    sed -e "s/PREFIX/sim-far-away-0fixed.${mode}/g" analysis-scripts/input.template.xml > input.xml
    sed -i -e "s/^.*FIXATOMS.*$//g" input.xml
    sed -i -e "s/MODE/${mode}/g" input.xml
    sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
    sed -i -e "s:init\.xyz:init.far-away.xyz:g" input.xml
    sed -i -e "s:in\.lmp:in.far-away.lmp:g" input.xml
    sed -i -e "s:data\.water:data.far-away.water:g" input.xml
    subdir="${mode}-far-away-0fixed"
    echo $subdir
    if [ ! -d $subdir ]; then mkdir $subdir; fi
#    cp {input.xml,data.far-away.water,in.far-away.lmp,init.far-away.xyz} $subdir/
    cp input.xml $subdir/

    # Far away, 1 fixed atom
    sed -e "s/PREFIX/sim-far-away-1fixed.${mode}/g" analysis-scripts/input.template.xml > input.xml
    sed -i -e "s/FIXATOMS/3/g" input.xml
    sed -i -e "s/MODE/${mode}/g" input.xml
    sed -i -e "s:LMP_EXE:${LMP_EXE}:g" input.xml
    sed -i -e "s:init\.xyz:init.far-away.xyz:g" input.xml
    sed -i -e "s:in\.lmp:in.far-away.lmp:g" input.xml
    sed -i -e "s:data\.water:data.far-away.water:g" input.xml
    subdir="${mode}-far-away-1fixed"
    echo $subdir
    if [ ! -d $subdir ]; then mkdir $subdir; fi
#    cp {input.xml,data.far-away.water,in.far-away.lmp,init.far-away.xyz} $subdir/
    cp input.xml $subdir/
done

for dir in *-?fixed; do cd $dir; rm data.water; ln -s ../data.water ./; rm init.xyz; ln -s ../init.xyz; rm in.lmp; ln -s ../in.lmp; cd ..; done
for dir in *far-away-?fixed; do cd $dir; rm data.water; rm init.xyz; rm in.lmp; cd ..; done
for dir in *far-away-?fixed; do cd $dir; rm data.far-away.water; ln -s ../data.far-away.water ./; rm init.far-away.xyz; ln -s ../init.far-away.xyz; rm in.far-away.lmp; ln -s ../in.far-away.lmp; cd ..; done

rm input.xml
