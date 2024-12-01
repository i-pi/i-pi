#Define model parameters
  eta=1.5

  wb=500
  wc=500
  V0=2085
  mass=1837.36223469
  epsilon=-1.0
  delta=0
  epsilon2=0
  deltaQ=1

  address=localhost
  model='DW_friction'


#Launch i-pi and wait
  i-pi input.xml & 
  sleep 3

#Launch driver
   #i-pi-driver-py -m ${model} -o ${wb},${V0},${mass},${x0},${eta},${epsilon},${delta},${deltaQ} -u -a ${address} 
   arg="w_b=${wb},v0=${V0},m=${mass},delta=${delta},eta0=${eta},eps1=${epsilon},eps2=${epsilon2},deltaQ=${deltaQ}"

   echo ${arg}
#   def __init__(self,w_b=None,v0=None,m=None,eta0=None,eps1=None,eps2=None,delta=None,deltaQ=None, *args, **kwargs):
   i-pi-driver-py -m ${model} -o ${arg} -u -a ${address} 
