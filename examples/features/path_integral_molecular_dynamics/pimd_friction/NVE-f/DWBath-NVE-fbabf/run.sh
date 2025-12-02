set -e -x

#Define model parameters
#This parameters are cpoied from the ring polymer instanton examples (/i-pi/examples/features/ring_polymer_instanton/rates_1D_double_well_with_friction/implicit_bath_pos-independent/run.sh). 
eta=1.5
wb=500
wc=500
V0=2085
mass=1837.36223469
x0=0
epsilon=0
delta=0
deltaQ=1

address=localhost
model='DW_bath'
ipi="i-pi"
#Launch i-pi and wait
$ipi input.xml > i-pi.out &
sleep 3

#Launch driver
i-pi-driver-py -m ${model} -o ${wb},${V0},${mass},${x0},${eta},${epsilon},${delta},${deltaQ},${wc} -u -a ${address} > DW-bath.out &

wait

echo "# Simulation complete"
