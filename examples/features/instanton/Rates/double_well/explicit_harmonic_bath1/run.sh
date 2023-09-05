
#Launch i-pi and wait
       i-pi input.xml &
       sleep 2
#Laun Driver
       #i-pi-driver -u -m harmonic_bath -o 1,1.674,500
       i-pi-driver -u -m harmonic_bath -o 1,0.5,500
  
