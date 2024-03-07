ipi=i-pi
driver=i-pi-driver
sleep_time=2

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing driver"
sleep ${sleep_time}

${driver} -u -a base -m qtip4pf &
${driver} -u -a comm-1 -m qtip4pf-c-1-delta &
${driver} -u -a comm-2 -m qtip4pf-c-2-delta &

echo "# All driver instances are running"

wait

echo "# Simulation complete"
