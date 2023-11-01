ipi=i-pi
driver=i-pi-driver 
sleep_time=10

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing driver"
sleep ${sleep_time}

${driver} -m qtip4pf -u -a qtip4pf > /dev/null & 
echo "# driver is running"

wait

echo "# Simulation complete"
