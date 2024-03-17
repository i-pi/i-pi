ipi=i-pi
sleep_time=10

cmd="${ipi} input.xml  > i-pi.out &"
echo ${cmd}
eval ${cmd}
echo "# i-PI is running"

echo "# Waiting for ${sleep_time} (s) before executing pswater"
sleep ${sleep_time}

cmd="i-pi-driver -a driver -u -m pswater  > /dev/null &"
echo ${cmd}
eval ${cmd}
echo "# pswater is running"

wait

echo "# Simulation complete"