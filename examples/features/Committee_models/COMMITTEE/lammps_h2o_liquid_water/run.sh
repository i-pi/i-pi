ipi=i-pi
lmp=lmp_mpi 
sleep_time=10

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing LAMMPS"
sleep ${sleep_time}

${lmp} < in-1.lmp > /dev/null & 
${lmp} < in-2.lmp > /dev/null & 
${lmp} < in-3.lmp > /dev/null & 
echo "# LAMMPS instances are running"

wait

echo "# Simulation complete"
