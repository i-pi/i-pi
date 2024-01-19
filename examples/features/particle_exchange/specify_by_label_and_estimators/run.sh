ipi=i-pi
driver="i-pi-driver -m harm3d -o 1.21647924E-8 -u -h bosons-trapped"
sleep_time=4

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing driver"
sleep ${sleep_time}

${driver} > /dev/null & 
echo "# Driver is running"

wait

echo "# Simulation complete"
