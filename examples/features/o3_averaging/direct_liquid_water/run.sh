ipi=i-pi
driver1="i-pi-driver -m qtip4pf -u -a h2o-base"
sleep_time=4

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing driver"
sleep ${sleep_time}

${driver1} > /dev/null & 
echo "# Driver is running"

wait

echo "# Simulation complete"
