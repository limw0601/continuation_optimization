orbit_transfer_7p corresponds to the case releasing the two ending points of transfer trajectory at the last
two stages

orbit_transfer_new_path corresponds to the case allowing the two ending points free to move in the first stage

in orbit_transfer_7p, I validated the code by solving the Hamiltonian two-point boundary-value problem directly
using continuation methods (see HTBVP, for the case with fixed ending points). The optimal state matches well
and the control has small difference at the initial time due to the 10-term expansion in our framework

in two cases, please run demo_auto_het_instances as the main program

a typo in paper, for the initial design, J=8.3237 INSTEAD OF J=8.3273