ipi=i-pi
driver="i-pi-driver -m pswater -u -a h2o-geop"
sleep_time=4

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing driver"
sleep ${sleep_time}

${driver} > /dev/null & 
echo "# Driver is running"

wait

echo "# Simulation complete"
