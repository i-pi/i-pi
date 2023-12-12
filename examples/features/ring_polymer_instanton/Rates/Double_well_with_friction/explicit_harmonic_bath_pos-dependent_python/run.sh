#Define model parameters
  eta=1.5

  wb=500
  wc=500
  V0=2085
  mass=1837.36223469
  x0=0
  epsilon=-1.0
  delta=0
  deltaQ=1

  address=localhost
  model='DW_bath'


#Launch i-pi and wait
  i-pi input.xml &
  sleep 3

#Launch driver
   i-pi-driver-py -m ${model} -o ${wb},${V0},${mass},${x0},${eta},${epsilon},${delta},${deltaQ},${wc} -u -a ${address}
