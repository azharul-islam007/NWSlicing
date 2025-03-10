function [] = test_slicing_RL_DL_f()
close all;
% initialization of random seed for repeatibility  
rand_seed = 300;  

% number of episodes. Every 'action' is evaluated in one episode and a new
% 'action' will be taken after the completition of an episode.
nr_episode = 100; 

% Temperature parameter, refer to section 'SELECTION CRITERION' for its
% impact on selection probability: with a high value of tao, te action
% probabilities become nearly equal and it will take longer time for the
% action selection in RL algorithm to converge; with a low value of tao,
% the action selection will put on more weights on the Q-value, and it will
% take shorter time for convergence.
tao = 1; 

% The averaging time windows, T, in Table 1, in unit of slot.
% In the paper it is recommended to use 3s/0.5e-3=6000, in which, though RL
% performance will be much robust, it will take much much longer time to
% complete the simulation. 
% In this simulation in a laptop, it is recommemend to use use smaller
% value, such 10, 100.
Avg_T = 10; 

% The number of vehicular UEs in the four cluster in the off-line network
% model. It should be tuned according to the configuration of online
% network model (which could be simulated or could be in a field network)
V_v2x = [1 2 3 4] * 2;

% The number of the eMBB UEs in the off-line network model. It should also
% be tuned according to the configuration of online network model (which
% could be simulated or could be in a field network)
M_eMBB = 12 * 2;

% The number of the mMTC UEs in the off-line network model. It should also
% be tuned according to the configuration of online network model (which
% could be simulated or could be in a field network)
M_mMTC = 4000;

slicing_RL_DL_f(rand_seed, nr_episode, tao, Avg_T, V_v2x, M_eMBB,M_mMTC);
end

