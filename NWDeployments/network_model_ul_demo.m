close all; 
M_eMBB = 6; % Number of eMBB UE randomly distributed in the cell
M_veh = 4; % Number of vehicular randomly distributed in the cell in the highway
M_mMTC = 4; % Number of mMTC UE randomly distributed in the cell
lambda_e = 1; % eMBB session generate rate (session/second), refer to figure 3 in paper
sim_time = 10; % number of simulation drop
nw = Network_Model_ul_c(M_eMBB, M_veh, M_mMTC, lambda_e, sim_time); % create the network model for dl simulation
nw.start_simulation(); % start to run the network model
