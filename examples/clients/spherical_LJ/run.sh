# source /home/stoccoel/codes/i-pi/env.sh
i-pi input.xml > i-pi.out & 
sleep 5
i-pi-driver -u -a qtip4pf -m qtip4pf
