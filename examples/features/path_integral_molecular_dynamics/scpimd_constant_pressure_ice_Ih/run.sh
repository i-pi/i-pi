ipi=i-pi
lmp=lmp_base
sleep_time=10

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing LAMMPS"
sleep ${sleep_time}

${lmp} < in.lmp > /dev/null & 
echo "# LAMMPS is running"

wait

echo "# Simulation complete"
