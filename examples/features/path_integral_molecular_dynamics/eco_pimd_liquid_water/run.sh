ipi=i-pi
driver='i-pi-driver -u -a h2o-pimd.4 -m qtip4pf'
sleep_time=4

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing DRIVER"
sleep ${sleep_time}

${driver} & 
echo "# DRIVER is running"

wait

echo "# Simulation complete"
