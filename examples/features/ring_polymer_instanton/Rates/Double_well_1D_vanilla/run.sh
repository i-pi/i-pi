#Define model parameters

  wb=500
  V0=2085
  mass=1837.36223469
  x0=0.00
  address=localhost
  model='DoubleWell'


#Launch i-pi and wait
  i-pi input.xml & 
  sleep 3

#Launch driver
   i-pi-driver-py -m ${model} -o ${wb},${V0},${mass},${x0} -u -a ${address} 
