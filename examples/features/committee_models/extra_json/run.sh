ipi=i-pi
driver=i-pi-driver
sleep_time=2

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing driver"
sleep ${sleep_time}

${driver} -u -a h2o-comm.json -m qtip4pf-c-json &

echo "# All driver instances are running"

wait

echo "# Simulation complete"
