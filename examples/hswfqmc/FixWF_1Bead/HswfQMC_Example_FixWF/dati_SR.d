&dati_SR
SR_kind="LagrDynam"             !kind of SR algorithm used to determine the change of variational parameters: "plain_SR_", "fast_SR__", "Hgrad_min", "LagrDynam", ["Lagr_molt" not working], ["LinMolLag" not working]
SR_num_max_WO_MIN=-1            !number of maximum SR steps without a new minimum
SR_beta=0.000001               !initial SR_beta
SR_beta_Rp=0.000001            !(deprecated) initial SR_beta_Rp (for protonic derivatives)
SR_maxdeltaPsi=0.1             !maximum delta psi allowed for a SR step (assuming psi is normalized to 1
SR_max_SVD_MIN=0.000001        !maximum Singula Value Decomposition minimum value to consider when inverting the SR matrix
SR_change_bound=T              !set some bounds on the change of the variational parameters ?
SR_min_change=0.               !minimum change (in percentage) of the variational parameters in a SR step
SR_max_change=10.              !maximum change (in percentage) of the variational parameters in a SR step
SR_change_bound_Rp=F           !set some bounds on the change of the variational parameters ?
SR_min_change_Rp=0.            !minimum change (in percentage) of the variational parameters in a SR step
SR_max_change_Rp=10.           !maximum change (in percentage) of the variational parameters in a SR step
SR_adaptative_beta=F           !use the adaptative beta scheme ? (increase beta when SR does not found a new minimum for too long)
SR_lambda=T                    !use the lambda adaptative scheme ? (explained in reference ...)
lambda_init=1.                 !initial value for lambda
min_lambda=0.001                !minimum value for lambda
max_lambda=100.                 !maximum value for lambda
SR_lambda_Rp=F                 !use the lambda adaptative scheme for the protonic coordinates ?
lambda_Rp_init=1.              !initial value for lambda_Rp
min_lambda_Rp=0.001             !minimum value for lambda_Rp
max_lambda_Rp=100.              !maximum value for lambda_Rp
/
