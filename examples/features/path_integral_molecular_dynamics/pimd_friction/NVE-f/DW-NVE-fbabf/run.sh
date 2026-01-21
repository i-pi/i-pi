set -e -x

#Define model parameters
eta=1.5
wb=500
wc=500
V0=2085
mass=1837.36223469
delta=0
eps1=0
eps2=0
deltaQ=1

address=localhost
model='DW_friction'
ipi="i-pi"
#Launch i-pi and wait
$ipi input.xml > i-pi.out &
sleep 3

#Launch driver
arg="w_b=${wb},v0=${V0},m=${mass},delta=${delta},eta0=${eta},eps1=${eps1},eps2=${eps2},deltaQ=${deltaQ}"
i-pi-driver-py -m ${model} -o ${arg} -u -a ${address} > DW.out &

wait

echo "# Simulation complete"
