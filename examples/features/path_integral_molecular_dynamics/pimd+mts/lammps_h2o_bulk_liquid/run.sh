ipi=i-pi
lmp=lmp_mpi 
sleep_time=10

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing LAMMPS"
sleep ${sleep_time}

${lmp} < in_short.lmp > /dev/null & 
${lmp} < in_long.lmp > /dev/null & 
echo "# LAMMPS is running"

wait

echo "# Simulation complete"
