ipi=i-pi
lmp=lmp_mpi 
sleep_time=10

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing LAMMPS"
sleep ${sleep_time}

${lmp} < in_baseline.lmp > /dev/null & 
${lmp} < in_delta_1.lmp > /dev/null & 
${lmp} < in_delta_2.lmp > /dev/null & 
echo "# LAMMPS instances are running"

wait

echo "# Simulation complete"
