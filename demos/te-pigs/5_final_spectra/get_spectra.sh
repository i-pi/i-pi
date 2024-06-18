IPI_GETACF=~/source/i-pi/bin/i-pi-getacf 

${IPI_GETACF} -ifile ../4_dielectric_response_prediction/positions.xyz -mlag 2048 -ftpad 2048 -ftwin cosine-hanning -oprefix xx -der 
