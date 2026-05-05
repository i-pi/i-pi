ipi=i-pi
sleep_time=4

${ipi} input_multiple_drivers.xml > log.i-pi & 
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing driver"
sleep ${sleep_time}

i-pi-driver -m qtip4pf -u -a h2o-cl-cavmd > /dev/null & 
i-pi-driver -u -h h2o-dipole-cl-cavmd -m water_dip_pol -o 0 &> out
echo "# Driver is running"

wait

echo "# Simulation complete"
