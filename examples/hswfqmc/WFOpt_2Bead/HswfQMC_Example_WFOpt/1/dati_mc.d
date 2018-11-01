&dati_mc
N_mc=-10000                                     !Number of sampled points (or steps). If N_mc > 0 then it specifies the number for each processor
N_blank=1000                                     !Number of warming steps
N_1ppt=-150                                      !Number of performed single-particle moves between two consecutive sampled points (or steps). If N_1ppt < 0, then N_1ppt=-N*N_1ppt/100 (Np is the number of particles)
flag_TABC=T                                      !Twist Averaged Boundary Conditions (T=true or F=false)
N_TABC=1000                                      !After how many sampled points a twist is performed. If N_TABC < 0, there will be -N_TABC twist/processor
N_mc_relax_TABC=100                              !Number of relaxing steps after a twist
step_e=0.5                                       !Initial step length for the electrons
step_se=0.1                                      !Initial step length for the electronic shadows
step_p=0.075                                     !(NOT IMPLEMENTED YET) Initial step length for the protons 
acceptance_rate=50                               !Target acceptance rate. The step lengths will be automaticall adjusted accordingly
flag_continua=F                                  !Resume a previous simulation (T or F)
howtomove='1ppt'                                 !Walking method: 'allp'=all together    '1ppt'=single-particle moves
propmove='gaus'                                  !Transition pdf: 'flat'=flat distribution   'gaus'=gaussian e^(-8*(x-x_0)^2/step^2)
trimer_steps=F                                   !Propose Trimers (electron-shadow1-shadow2) moves. Important if the Kernel is very tight 
flag_elettroni=T                                 !Integrate over the electronic coordinates (T or F)
flag_protoni=F                                   !(NOT IMPLEMENTED YET) Integrate over the protonic coordinates
flag_shadow=F                                    !Integrate over the electronic-shadow coordinates
flag_E_tot=T                                     !Compute the totale energy
flag_E_kin=F                                     !Compute the kinetic energy (Pandharipande-Bethe and Jackson-Feenber)
flag_E_pot=F                                     !Compute the potential energy
flag_somme_ewald=T                               !Use the Ewald summation
alpha_ewald=-1.d0                               !alpha parameters used in the Ewald summation. If it's equal to -1.d0, then it is set automatically
num_k_ewald=515                                  !Number of k vectors used in the Ewald summation
flag_gr=F                                        !Compute the pair correlation functions
N_hist=250                                       !Number of points in the grid used for the pair correlation functions
flag_posizioni=F                                 !Save into files located in positions/<...>.pos the final positions of the walkers (necessary for resuming the calculation later on) 
flag_normalizza_pos=T                            !Normalize the position coordinates so that they stay in the range [0,1]
N_AV=-6400                                       !Number of data averaged before being written into file. If N_AV < 0, then data are packed into -N_AV elements (a good value is -6400)
flag_MPI=T                                       !Use MPI
what_to_do='stocrec'                             !Which task has to be performed: 'simpcal'=VMC calculation   'stocrec'=stochastic reconfiguration minimization (others are deprecated)
stampa_dati_funzione_onda=T                      !Print informations about the wave function during minimization
path_dati_funzione_onda='wf_now.d'               !Path that specifies the file containing the informations about the trial wave function to be used
accuracy_energy_opt=0.001d0                      !(DEPRECATED)
flag_disk=F                                      !Save all the data on files for computing the mean and variance. If not, everything is made only using the RAM
flag_output=T                                    !Save all the outputs on the output.d file
quick_error=16		                         !Number of blocks used to estimate the error with the blocking technique. If quick_error=0 then the optimal value is determined automatically
flag_random_file=T                               !Use a file with true random numbers to initialize the pseudo-random number generator
random_seed_path='randomseed.d'                 !Path to the file used to initialized the pseudo-random number generator
/


