ipi=i-pi
cmd="${ipi} input.xml  > i-pi.out &"
echo ${cmd}
eval ${cmd}
echo "# i-PI is running"
wait
echo "# Simulation complete"