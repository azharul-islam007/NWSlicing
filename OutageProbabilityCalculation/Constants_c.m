classdef Constants_c < handle
    properties (Constant)
        %% The following parameters could be tuned for performance comparison
        % For performance comparision when the offline model is NOT
        % utilizing the online model's inputs:
        ignore_online_nw_characteristics = false;

        % Number of vehicular randomly distributed in the cell in the highway
        M_veh = 8; 

        % V2X packet arrval rate (packets/second)
        lambda_niu = 200;    

        % Number of eMBB UE randomly distributed in the cell
        M_eMBB = 4; 

        % lambda_e: eMBB session generate rate (session/second), refer to figure 3 in paper
        lambda_e_lst = 300 : 100 : 600;            

        % Number of mMTC UE randomly distributed in the cell
        M_mMTC = 12; 

        % the wakeup-transmission probability of an mMTC packet
        p_mMTC_tx = 0.75;  

        % the max packet size of mMTC packets (bytes)
        max_mMTC_packet_size = 128;        
        
        % number of simulation drop in each round
        sim_nr_drops_per_round = 100; 

        % number of simulation round
        sim_nr_round = 20; 

        % eMBB UE arrival rate, 1UE/second
        lambda_m = 1; 

        % V2X vehicular arrival rate, 1 UE/second
        lambda_a = 1; 
        
        % Number of action for slicing of mMTC
        Ar = 4; 

        % Number of actions for slicing allocation between eMBB and V2X
        Ax = 20; 

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
        % In this simulation in a laptop, it is recommemend to use a smaller
        % value, such 10, 100.
        Avg_T = 20;    

        % RL algorithm learning rate
        alpha_learning_rate = 0.1;

        % 
        % Use a normal distribution for SINR dB by default
        pd_SINR_dB_by_default = makedist('Normal', 'mu', 20, 'sigma', 20);

        % Use a Poission distribution for number of packets per drop for V2X, mMTC traffic by default
        pd_Nr_packet_by_default = makedist('Poisson', 'lambda', 1);   

        % Use a Poission distribution for number of session per drop for eMBB traffic by default
        pd_Nr_session_by_default = makedist('Poisson', 'lambda', 1);       

        % Use a binomial distribution for the transmission probability for
        % mMTC traffic by default
        pd_tx_by_default = makedist('Binomial', 'N', 1, 'p', 0.5);    

        % Use a binomial distribution for the transmission probability for
        % mMTC traffic by default
        pd_log2_packet_size_by_default = makedist('Multinomial', 'Probabilities', ones(7, 1) / 7);         

        %% The following parameters SHOULD NOT be tunned
        % packet size for V2X service
        Sm = 300; 

        % velocity in meter/second
        velocity = 80/3.6; 

        % probability of selection between sidelink or cellular in the UL
        % for a vehicular UE
        p_SL = 0.5; 

        % number of PRB, which is the same for UL and DL
        N_RB = 200; 

        % Subcarrier bandwidth
        f_sc = 30e3; % 30KHz

        % one PRB consists of 12 subcarriers
        N_sc_RB = 12; 

        % 4 clusters uniformally distributed along the highway
        nr_clusters = 4; 

        % 0.5ms, the time duration of a simulation drop
        T_drop = 0.5e-3;  

        % 1Mbps, guaranteed bit rates for eMBB sessions
        Rb = 1e6; 

        % TTI duration, should be 0.5ms because \Delta_f (subcarrier interval) is 30kHz
        Fd = 0.5e-3; 

        % 3 service types, {'V2X', 'eMBB', 'mMTC'}
        Nr_S = 3; 

        % The first service type is V2X
        I_V2X = 1; 

        % The second service type is eMBB
        I_eMBB = 2; 

        % The third service type is mMTC
        I_mMTC = 3;         

        % one eMBB session requires 1Mbps throughput
        Rb_session = 1e6; 

        % bandwidth of a PRB, equal to 30kHz * 12 subcarriers
        B = 30e3 * 12; 

        % cell radius
        radius = 500; 

        % cell diameter
        diameter = 1000; 

        % the height of the antenna at the gNB base station
        gNB_Bs_height = 10; 

        % the height of the antenna at the UE or the car
        UE_height = 1.5; 

        % the height of the rsu
        rsu_height = 3; 
        
        %% The following parameters are RECOMMENDED and should be kept as it is

        % a distance between the gNB and the highway
        d_gNB_highway = 30; 

        % a distance between the Rsu and the highway
        d_Rsu_highway = 25; 

        % the width of a lane in the highway
        d_lane = 24; 

        % number of lanes in one direction
        nr_lane_one_direction = 3; 

        % number of lanes in bi directions, it should be equal to the double of that in one direction
        nr_lane_bi_direction = Constants_c.nr_lane_one_direction * 2; 

        % position x: the gNB Bs should be at the center of the hexagonal region
        gNB_Bs_pos_x = 500; 

        % position y: the gNB Bs should be at the center of the hexagonal region
        gNB_Bs_pos_y = 433 + 20 + 3 * 24;  

        % eMBB UEs are distributed at a distance not less than radius_min
        radius_min = 25; 

        % eMBB UEs are distributed at a distance not larger than radius_max
        radius_max = Constants_c.radius / 2 * sqrt(3); 

        % eMBB session minimum duration time
        T_session_eMBB_min = 100; 

        % eMBB session maximum duration time
        T_session_eMBB_max = 140; 

        % The probability of selecting sidelink against cellular mode.
        % (Default: equally)
        p_sidelink_cellular_selection = 0.5;     

        % log file Id to save the matlab console globally
        fid = fopen('matlab_console.log', 'w+');

    end

    methods(Static)
        function [] = print_configurations()
            fprintf(Constants_c.fid, '---------------Configurations:------------------ \n');
            fprintf(Constants_c.fid, 'ignore_online_nw_characteristics = %s\n', string(Constants_c.ignore_online_nw_characteristics));
            fprintf(Constants_c.fid, 'M_veh = %d\n', Constants_c.M_veh); 
            fprintf(Constants_c.fid, 'lambda_niu = %f\n', Constants_c.lambda_niu);    
            fprintf(Constants_c.fid, 'M_eMBB = %d\n', Constants_c.M_eMBB);         
            fprintf(Constants_c.fid, 'M_mMTC = %d\n', Constants_c.M_mMTC); 
            fprintf(Constants_c.fid, 'p_mMTC_tx = %f\n', Constants_c.p_mMTC_tx);  
            fprintf(Constants_c.fid, 'max_mMTC_packet_size = %d\n', Constants_c.max_mMTC_packet_size);        
            fprintf(Constants_c.fid, 'sim_nr_drops_per_round per round = %d\n', Constants_c.sim_nr_drops_per_round); 
            fprintf(Constants_c.fid, 'sim_nr_round = %d\n', Constants_c.sim_nr_round); 
            fprintf(Constants_c.fid, '------------------------------------------------ \n');
        end

		% ! Remember to close the globally open fid
        function [] = close_fid()
            fclose(Constants_c.fid);
        end
    end

end