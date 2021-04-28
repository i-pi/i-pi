n=${1}
tail -n +2 inst.instanton_FINAL_${n}.ener |awk '{print $2}' > aux1
grep H inst.instanton_FINAL_${n}.xyz |awk '{print $2}' >aux2
paste aux2 aux1 >inst1D.dat
rm aux1 aux2
xmgrace  inst1D.dat ../REF/friction_data/MEP.dat
