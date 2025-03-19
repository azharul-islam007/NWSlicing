classdef Constants_c < handle
    properties (Constant)
        radius = 500; % cell radius
        diameter = 1000; % cell diameter
        gNB_Bs_height = 10; % the height of the antenna at the gNB base station
        UE_height = 1.5; % the height of the antenna at the UE or the car
        rsu_height = 3; % the height of the rsu
        d_gNB_highway = -48; % a distance between the gNB and the highway
        d_Rsu_highway = 24; % a distance between the Rsu and the highway
        d_lane = 16; % the width of a lane in the highway
        nr_lane_one_direction = 3; % number of lanes in one direction
        nr_lane_bi_direction = 6; % number of lanes in bi directions, it should be equal to the double of that in one direction
        gNB_Bs_pos_x = 500; % position x: the gNB Bs should be at the center of the hexagonal region
        gNB_Bs_pos_y = 433;  % position y: the gNB Bs should be at the center of the hexagonal region
        T_drop = 0.5e-3; % 0.5ms, the time duration of a simulation drop 
        Rb = 1e6; % 1Mbps, guaranteed bit rates for eMBB sessions
        T_session_eMBB_min = 100; % session duration time (minimum)
        T_session_eMBB_max = 140; % session duration time (maximum)
        Sm = 300; % packet size for V2X service
        velocity = 80/3.6; % kmh -> meter/second
        p_SL = 0.5; % probability of selection between sidelink or cellular in the UL for a vehicular UE
        N_RB = 200; % number of PRB
        f_sc = 30e3; % 30KHz
        N_sc_RB = 12; % one PRB consists of 12 subcarriers
        nr_clusters = 4; % 4 clusters uniformally distributed along the highway
    end

end