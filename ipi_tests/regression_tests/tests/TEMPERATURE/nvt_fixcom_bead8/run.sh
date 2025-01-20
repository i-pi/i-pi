ipi=i-pi
driver='i-pi-driver -u -a h2o-pimd -m qtip4pf'
sleep_time=4

${ipi} input_pileg.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing DRIVER"
sleep ${sleep_time}

${driver} &
${driver} &
${driver} &
${driver} & 
echo "# DRIVER is running"

wait

${ipi} input_pilel.xml > log.i-pi &
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing DRIVER"
sleep ${sleep_time}

${driver} &
${driver} &
${driver} &
${driver} &
echo "# DRIVER is running"

echo "# Simulation complete"
