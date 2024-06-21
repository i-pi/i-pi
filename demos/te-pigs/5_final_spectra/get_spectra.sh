IPI_GETACF=~/source/i-pi/bin/i-pi-getacf

# VDOS
${IPI_GETACF} -ifile ../3_production_simulations/simulation.pos_0.extxyz -mlag 2048 -ftpad 2048 -ftwin cosine-hanning -oprefix xx -dt "2.0 femtosecond"  -der

# IR
awk '{print 1; print "#"; print "H", $0}' ../4_dielectric_response_prediction/MACE_mu.txt > tmp.xyz
${IPI_GETACF} -ifile tmp.xyz -mlag 2048 -ftpad 2048 -ftwin cosine-hanning -oprefix mm -dt "2.0 femtosecond" -der
rm tmp.xyz

# Isotropic Raman
awk '{print 1; print "#"; print "H", $1 / 3, $5 / 3, $9 / 3}' ../4_dielectric_response_prediction/MACE_alpha.txt > tmp.xyz
${IPI_GETACF} -ifile tmp.xyz -mlag 2048 -ftpad 2048 -ftwin cosine-hanning -oprefix L0L0 -dt "2.0 femtosecond"  -der
rm tmp.xyz

# Anisotropic Raman
awk '{print 2; print "#"; print "H", $2, $3, $4; print "H", $5, $6, $0}' ../4_dielectric_response_prediction/MACE_alpha_sh.txt > tmp.xyz
${IPI_GETACF} -ifile tmp.xyz -mlag 2048 -ftpad 2048 -ftwin cosine-hanning -oprefix L2L2 -dt "2.0 femtosecond" -der
rm tmp.xyz
