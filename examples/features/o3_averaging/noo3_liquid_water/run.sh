ipi=i-pi
driver1="i-pi-driver -m noo3-h2o -u -a h2o-noo3"
driver2="i-pi-driver -m qtip4pf -u -a h2o-base"
sleep_time=4

${ipi} input.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing driver"
sleep ${sleep_time}

${driver1} > /dev/null & 
${driver2} > /dev/null &
echo "# Driver is running"

wait

echo "# Simulation complete"
