classdef Network_Model_ul_c < Network_Model_c
    properties
        Rsu_clusters;  % Four RSUs are deployed in every one of the four corners in every clusters
        M_rsu; % number of RSUs, should be equal to: num_of_corners_per_cluster (4) * num_of_cluster (4)
        veh_UEs_in_cellular; % the list of the vehicular UEs using cellular mode for UL transmission
        veh_UEs_in_sidelink; % the list of the vehicular UEs using sidelink mode for UL multicast within the cluster
        pd_veh_celluar_or_sidelink; % a distribution to simulate the UE's decision on using sidelink or cellular mode

        M_veh_UE_cellular; % number of UE working in cellular mode
        M_veh_UE_sidelink; % number of UE working in sidelink mode

        nr_corners_per_cluster;

        sim_nr_round;
        avg_prb_utilization_ratio_v2x; % average PRB utilization separately for V2X, eMBB, mMTC services in each round
        avg_prb_utilization_ratio_eMBB;
        avg_prb_utilization_ratio_mMTC;
        avg_prb_utilization_ratio;     

        nw_characteristics_v2x_CL;
        nw_characteristics_v2x_SL;

        avg_outage_rate_v2x; % average outage rate separately for V2X, eMBB, mMTC services in each round
        avg_outage_rate_eMBB;
        avg_outage_rate_mMTC;
        avg_outage_rate;
    end

    methods
        function self = Network_Model_ul_c(M_eMBB_i, M_veh_i, M_mMTC_i, lambda_e_i, sim_nr_drops_per_round_i, sim_nr_round_i)
            self@Network_Model_c(M_eMBB_i, M_veh_i, M_mMTC_i, lambda_e_i, sim_nr_drops_per_round_i);

            % Assuming UE can select either sidelink or cellular mode equally
            self.pd_veh_celluar_or_sidelink = makedist("Binomial", "N", 1, "p", Constants_c.p_sidelink_cellular_selection); 

            self.sim_nr_round = sim_nr_round_i;
            self.avg_prb_utilization_ratio = zeros(self.sim_nr_round, 1);
            self.avg_prb_utilization_ratio_v2x = zeros(self.sim_nr_drops_per_round, 1);
            self.avg_prb_utilization_ratio_eMBB  = zeros(self.sim_nr_drops_per_round, 1);
            self.avg_prb_utilization_ratio_mMTC = zeros(self.sim_nr_drops_per_round, 1);

            self.avg_outage_rate_v2x = zeros(self.sim_nr_round, 1);
            self.avg_outage_rate_eMBB = zeros(self.sim_nr_round, 1);
            self.avg_outage_rate_mMTC = zeros(self.sim_nr_round, 1);
            self.avg_outage_rate = zeros(self.sim_nr_round, 1);
            
        end

        % @Override
        function [] = collect_nw_characteristic(self, t_i)  
            % Override the function definition in its subclass, becasue of
            % the special handling of cellular mode or sidelink mode
            % vehicular UE
            for count = 1 : self.M_veh_UE_cellular
                self.nw_characteristics_v2x_CL.SINRdB{t_i}{count} = self.veh_UEs_in_cellular{count}.sinr_dB;
                self.nw_characteristics_v2x_CL.nr_packets{t_i}{count} = self.veh_UEs_in_cellular{count}.nr_packet;
            end
            for count = 1 : self.M_veh_UE_sidelink
                self.nw_characteristics_v2x_SL.SINRdB{t_i}{count} = self.veh_UEs_in_sidelink{count}.sinr_dB;
                self.nw_characteristics_v2x_SL.nr_packets{t_i}{count} = self.veh_UEs_in_sidelink{count}.nr_packet;
            end            
            for count = 1 : self.M_eMBB
                self.nw_characteristics_eMBB.SINRdB{t_i}{count} = self.eMBB_UEs{count}.sinr_dB;
                self.nw_characteristics_eMBB.nr_session{t_i}{count} = self.eMBB_UEs{count}.nr_session();
            end
            for count = 1 : self.M_mMTC
                self.nw_characteristics_mMTC.SINRdB{t_i}{count} = self.mMTC_UEs{count}.sinr_dB;
                self.nw_characteristics_mMTC.packet_size{t_i}{count} = self.mMTC_UEs{count}.packet_size;
                self.nw_characteristics_mMTC.tx_flags{t_i}{count} = self.mMTC_UEs{count}.tx_flag;
            end
        end        

        function [] = init_simulation(self)
		% Beside the initialization in the basic network model, the RSUs at the corners of the clusters should also be created
            self.init_simulation@Network_Model_c();
            self.Rsu_clusters = cell(Constants_c.nr_clusters, 1);
            delta_x = Constants_c.diameter / (2 * Constants_c.nr_clusters);
            d_Rsu = Constants_c.d_Rsu_highway;
            delta_y = Constants_c.d_lane * Constants_c.nr_lane_one_direction;
            cluster_center_x = delta_x * (1 : 2 : 2 * Constants_c.nr_clusters);
            cluster_center_y = Constants_c.gNB_Bs_pos_y - Constants_c.d_gNB_highway - delta_y;
            for count = 1 : Constants_c.nr_clusters
                z = Constants_c.rsu_height;
                Rsu_nw = Rsu_c(self.next_ueId, cluster_center_x(count) - delta_x, cluster_center_y + delta_y + d_Rsu, z);
                Rsu_ne = Rsu_c(self.next_ueId + 1, cluster_center_x(count) + delta_x, cluster_center_y + delta_y + d_Rsu, z);
                Rsu_se = Rsu_c(self.next_ueId + 2, cluster_center_x(count) + delta_x, cluster_center_y - delta_y - d_Rsu, z);
                Rsu_sw = Rsu_c(self.next_ueId + 3, cluster_center_x(count) - delta_x, cluster_center_y - delta_y - d_Rsu, z);
                self.next_ueId = self.next_ueId + 4;
                self.Rsu_clusters{count} = [Rsu_nw; Rsu_ne; Rsu_se; Rsu_sw];
            end
            self.nr_corners_per_cluster = 4; % each sidelink UE is connected with 4 RSUs in the every corner.
            self.M_rsu = Constants_c.nr_clusters * self.nr_corners_per_cluster; 
        end

        function [] = prepare_simulation_round(self, r_i)
            self.nw_characteristics_v2x_CL = struct();
            self.nw_characteristics_v2x_CL.SINRdB = cell(self.sim_nr_drops_per_round, 1);
            self.nw_characteristics_v2x_CL.nr_packets = cell(self.sim_nr_drops_per_round, 1);

            self.nw_characteristics_v2x_SL = struct();
            self.nw_characteristics_v2x_SL.SINRdB = cell(self.sim_nr_drops_per_round, 1);
            self.nw_characteristics_v2x_SL.nr_packets = cell(self.sim_nr_drops_per_round, 1);

            self.nw_characteristics_eMBB = struct();
            self.nw_characteristics_eMBB.SINRdB = cell(self.sim_nr_drops_per_round, 1);
            self.nw_characteristics_eMBB.nr_session = cell(self.sim_nr_drops_per_round, 1);

            self.nw_characteristics_mMTC = struct();
            self.nw_characteristics_mMTC.SINRdB = cell(self.sim_nr_drops_per_round, 1);
            self.nw_characteristics_mMTC.packet_size = cell(self.sim_nr_drops_per_round, 1);
            self.nw_characteristics_mMTC.tx_flags = cell(self.sim_nr_drops_per_round, 1);

            self.nw_outage_v2x = zeros(self.sim_nr_drops_per_round, 1);
            self.nw_outage_eMBB = zeros(self.sim_nr_drops_per_round, 1);
            self.nw_outage_mMTC = zeros(self.sim_nr_drops_per_round, 1);

            self.prb_utilization_ratio_v2x = zeros(self.sim_nr_drops_per_round, 1);
            self.prb_utilization_ratio_eMBB = zeros(self.sim_nr_drops_per_round, 1);
            self.prb_utilization_ratio_mMTC = zeros(self.sim_nr_drops_per_round, 1);
            self.prb_utilization_ratio = zeros(self.sim_nr_drops_per_round, 1);
        end


        function [] = update_slicing_ratio(self, round_i)
            if (Constants_c.ignore_online_nw_characteristics)
                pd_v2x_CL_SINR_dB = Constants_c.pd_SINR_dB_by_default;
                pd_v2x_CL_Nr_packet = Constants_c.pd_Nr_packet_by_default;
                pd_v2x_SL_SINR_dB = Constants_c.pd_SINR_dB_by_default;
                pd_v2x_SL_Nr_packet = Constants_c.pd_Nr_packet_by_default;
                pd_eMBB_SINR_dB = Constants_c.pd_SINR_dB_by_default;
                pd_eMBB_nr_session = Constants_c.pd_Nr_session_by_default;
                pd_mMTC_SINR_dB = Constants_c.pd_SINR_dB_by_default;
                pd_mMTC_log2_packet_size = Constants_c.pd_log2_packet_size_by_default;
                pd_mMTC_tx = Constants_c.pd_tx_by_default;
            else
                [pd_v2x_CL_SINR_dB, pd_v2x_CL_Nr_packet, ...
                    pd_v2x_SL_SINR_dB, pd_v2x_SL_Nr_packet, ...
                    pd_eMBB_SINR_dB, pd_eMBB_nr_session, ...
                    pd_mMTC_SINR_dB, pd_mMTC_log2_packet_size, pd_mMTC_tx] = self.estimate_nw_parameters_f();
            end

            alpha = Constants_c.alpha_learning_rate; % RL learning rate
            RL_heuristic_ul = RL_heuristic_ul_c(alpha, pd_v2x_CL_SINR_dB, pd_v2x_CL_Nr_packet, ...
                pd_v2x_SL_SINR_dB, pd_v2x_SL_Nr_packet, ...
                pd_eMBB_SINR_dB, pd_eMBB_nr_session, ...
                pd_mMTC_SINR_dB, pd_mMTC_log2_packet_size, pd_mMTC_tx);

            rand_seed = round_i;
            M_veh_CL = ceil(self.M_veh * (1 - Constants_c.p_sidelink_cellular_selection));
            M_veh_SL = self.M_veh - M_veh_CL;
            V_v2x_CL = ones(1, Constants_c.nr_clusters) * floor(M_veh_CL / Constants_c.nr_clusters); % Assuming that all Veh UE are uniformally distributed in each cluster
            rem = mod(M_veh_CL, Constants_c.nr_clusters);
            % distribute the remaining UEs randomly in any one of the
            % clusters
            if (rem > 0)
                for i = 1 : rem
                    cluster_index = randsample(Constants_c.nr_clusters, 1);
                    V_v2x_CL(cluster_index) = V_v2x_CL(cluster_index) + 1;
                end
            end
            V_v2x_SL = ones(1, Constants_c.nr_clusters) * floor(M_veh_SL /  Constants_c.nr_clusters); % Assuming that all Veh UE are uniformally distributed in each cluster
            rem = mod(M_veh_SL, Constants_c.nr_clusters);
            % distribute the remaining UEs randomly in any one of the
            % clusters
            if (rem > 0)
                for i = 1 : rem
                    cluster_index = randsample(Constants_c.nr_clusters, 1);
                    V_v2x_SL(cluster_index) = V_v2x_SL(cluster_index) + 1;
                end
            end

            M_eMBB = self.M_eMBB;
            M_mMTC = self.M_mMTC;
            nr_episode = Constants_c.nr_episode;
            tao = Constants_c.tao;
            Avg_T = Constants_c.Avg_T;
            alpha_sel = RL_heuristic_ul.offline_RL_heuristic_learning_slicing_ratio(rand_seed, ...
                nr_episode, tao, Avg_T, V_v2x_CL, V_v2x_SL, M_eMBB, M_mMTC);            
            V2X = Constants_c.I_V2X;
            eMBB = Constants_c.I_eMBB;
            mMTC = Constants_c.I_mMTC;
            self.alpha_v2x = alpha_sel(V2X);
            self.alpha_eMBB = alpha_sel(eMBB);
            self.alpha_mMTC = alpha_sel(mMTC);
            fprintf(Constants_c.fid, 'slicing ratio updated after round (%d) = %.3f %.3f %.3f\n', round_i, self.alpha_v2x, self.alpha_eMBB, self.alpha_mMTC);
        end

        function [pd_v2x_CL_SINR_dB_o, pd_v2x_CL_Nr_packet_o, ...
                pd_v2x_SL_SINR_dB_o, pd_v2x_SL_Nr_packet_o, ...
                pd_eMBB_SINR_dB_o, pd_eMBB_nr_session_o, ...
                pd_mMTC_SINR_dB_o, pd_mMTC_log2_packet_size_o, pd_mMTC_tx_o] = estimate_nw_parameters_f(self)
            if (self.M_veh_UE_cellular > 1) % To fit the normal distribution, at least two elements are provided
                sinr_dB = cell2vec_f(self.nw_characteristics_v2x_CL.SINRdB);
                pd_v2x_CL_SINR_dB_o = fitdist(sinr_dB, 'Normal'); 
                nr_packets = cell2vec_f(self.nw_characteristics_v2x_CL.nr_packets);            
                pd_v2x_CL_Nr_packet_o = fitdist(nr_packets, 'Poisson');
            else
                pd_v2x_CL_SINR_dB_o = Constants_c.pd_SINR_dB_by_default;          
                pd_v2x_CL_Nr_packet_o = Constants_c.pd_Nr_packet_by_default;
            end

             if (self.M_veh_UE_sidelink > 1) % To fit the normal distribution, at least two elements are provided
                sinr_dB = cell2vec_f(self.nw_characteristics_v2x_SL.SINRdB);
                pd_v2x_SL_SINR_dB_o = fitdist(sinr_dB, 'Normal'); 
                nr_packets = cell2vec_f(self.nw_characteristics_v2x_SL.nr_packets);            
                pd_v2x_SL_Nr_packet_o = fitdist(nr_packets, 'Poisson');
             else
                pd_v2x_SL_SINR_dB_o = Constants_c.pd_SINR_dB_by_default;          
                pd_v2x_SL_Nr_packet_o = Constants_c.pd_Nr_packet_by_default;
             end

            sinr_dB = cell2vec_f(self.nw_characteristics_eMBB.SINRdB);
            pd_eMBB_SINR_dB_o = fitdist(sinr_dB, 'Normal'); 
            nr_session = cell2vec_f(self.nw_characteristics_eMBB.nr_session);            
            pd_eMBB_nr_session_o = fitdist(nr_session, 'Poisson');

            sinr_dB = cell2vec_f(self.nw_characteristics_mMTC.SINRdB);
            pd_mMTC_SINR_dB_o = fitdist(sinr_dB, 'Normal'); 

            packet_size = cell2vec_f(self.nw_characteristics_mMTC.packet_size);
            [GC, GR] = groupcounts(packet_size);
            log2_GR = log2(GR);
            max_log2_GR = max(log2_GR);
            prob = zeros(max_log2_GR, 1);
            for count = 1 : length(log2_GR)
                k = log2_GR(count);
                prob(k) = GC(count) / sum(GC);
            end
            pd_mMTC_log2_packet_size_o = makedist('Multinomial', 'Probabilities', prob);

            tx = cell2vec_f(self.nw_characteristics_mMTC.tx_flags);
            pd_mMTC_tx_o = fitdist(tx, 'Binomial', 'N', 1);           
        end        

        function [] = plot(self)
            [h_gNB_Bs_o, h_eMBB_UE_o, h_mMTC_UE_o] = self.plot@Network_Model_c();      

            % draw vehicular UEs sending UL in cellular mode
            nr_veh_UE_l2r = 0;
            nr_veh_UE_r2l = 0;
            for count = 1 : self.M_veh_UE_cellular
                if strcmpi(self.veh_UEs_in_cellular{count}.direction, 'l2r')
                    nr_veh_UE_l2r = nr_veh_UE_l2r + 1;
                elseif strcmpi(self.veh_UEs_in_cellular{count}.direction, 'r2l')
                    nr_veh_UE_r2l = nr_veh_UE_r2l + 1;
                end
            end
            veh_UE_pos_x = zeros(nr_veh_UE_l2r, 1);
            veh_UE_pos_y = zeros(nr_veh_UE_l2r, 1);
            if (nr_veh_UE_l2r > 0)

                count2 = 1;
                for count = 1 : self.M_veh_UE_cellular
                    if strcmpi(self.veh_UEs_in_cellular{count}.direction, 'l2r')
                        veh_UE_pos_x(count2) = self.veh_UEs_in_cellular{count}.pos_x;
                        veh_UE_pos_y(count2) = self.veh_UEs_in_cellular{count}.pos_y;
                        count2 = count2 + 1;
                    end
                end                
            end
            h_veh_UE_cellular_l2r_o = plot(veh_UE_pos_x, veh_UE_pos_y, ">k");

            veh_UE_pos_x = zeros(nr_veh_UE_r2l, 1);
            veh_UE_pos_y = zeros(nr_veh_UE_r2l, 1);
            if (nr_veh_UE_r2l > 0)

                count2 = 1;
                for count = 1 : self.M_veh_UE_cellular
                    if strcmpi(self.veh_UEs_in_cellular{count}.direction, 'r2l')
                        veh_UE_pos_x(count2) = self.veh_UEs_in_cellular{count}.pos_x;
                        veh_UE_pos_y(count2) = self.veh_UEs_in_cellular{count}.pos_y;
                        count2 = count2 + 1;
                    end
                end                
            end
            h_veh_UE_cellular_r2l_o = plot(veh_UE_pos_x, veh_UE_pos_y, "<m");

            for count = 1 : self.M_veh_UE_cellular
                x = [self.gNB_Bs.pos_x self.veh_UEs_in_cellular{count}.pos_x];
                y = [self.gNB_Bs.pos_y self.veh_UEs_in_cellular{count}.pos_y];
                plot(x, y, "--c");  
                text(x(2),y(2), sprintf('sinr=%.1f dB', self.vehicular_UEs{count}.sinr_dB));
            end

            % draw vehicular UEs sending UL in sidelink mode
            nr_veh_UE_l2r = 0;
            nr_veh_UE_r2l = 0;
            for count = 1 : self.M_veh_UE_sidelink
                if strcmpi(self.veh_UEs_in_sidelink{count}.direction, 'l2r')
                    nr_veh_UE_l2r = nr_veh_UE_l2r + 1;
                elseif strcmpi(self.veh_UEs_in_sidelink{count}.direction, 'r2l')
                    nr_veh_UE_r2l = nr_veh_UE_r2l + 1;
                end
            end
            veh_UE_pos_x = zeros(nr_veh_UE_l2r, 1);
            veh_UE_pos_y = zeros(nr_veh_UE_l2r, 1);
            if (nr_veh_UE_l2r > 0)
                count2 = 1;
                for count = 1 : self.M_veh_UE_sidelink
                    if strcmpi(self.veh_UEs_in_sidelink{count}.direction, 'l2r')
                        veh_UE_pos_x(count2) = self.veh_UEs_in_sidelink{count}.pos_x;
                        veh_UE_pos_y(count2) = self.veh_UEs_in_sidelink{count}.pos_y;
                        count2 = count2 + 1;
                    end
                end          
            end
            h_veh_UE_sidelink_l2r_o = plot(veh_UE_pos_x, veh_UE_pos_y, ">k");

            veh_UE_pos_x = zeros(nr_veh_UE_r2l, 1);
            veh_UE_pos_y = zeros(nr_veh_UE_r2l, 1);
            if (nr_veh_UE_r2l > 0)
                count2 = 1;
                for count = 1 : self.M_veh_UE_sidelink
                    if strcmpi(self.veh_UEs_in_sidelink{count}.direction, 'r2l')
                        veh_UE_pos_x(count2) = self.veh_UEs_in_sidelink{count}.pos_x;
                        veh_UE_pos_y(count2) = self.veh_UEs_in_sidelink{count}.pos_y;
                        count2 = count2 + 1;
                    end
                end
            end
            h_veh_UE_sidelink_r2l_o = plot(veh_UE_pos_x, veh_UE_pos_y, "<m");
                      
            for count = 1 : self.M_veh_UE_sidelink
                cluster = self.veh_UEs_in_sidelink{count}.cluster;
                Rsus = self.Rsu_clusters{cluster};
                for ri = 1 : self.nr_corners_per_cluster
                    Rsu = Rsus(ri);         
                    x = [Rsu.pos_x self.veh_UEs_in_sidelink{count}.pos_x];
                    y = [Rsu.pos_y self.veh_UEs_in_sidelink{count}.pos_y];
                    plot(x, y, "--r");  
                end
                text(x(2),y(2), sprintf('sinr=%.1f dB', self.veh_UEs_in_sidelink{count}.sinr_dB));
            end

            Rsu_pos_x = zeros(Constants_c.nr_clusters * self.nr_corners_per_cluster, 1);
            Rsu_pos_y = zeros(Constants_c.nr_clusters * self.nr_corners_per_cluster, 1);
            count2 = 1;
            for cluster = 1 : Constants_c.nr_clusters
                Rsus = self.Rsu_clusters{cluster};
                for ri = 1 : self.nr_corners_per_cluster
                    Rsu = Rsus(ri);         
                    Rsu_pos_x(count2) = Rsu.pos_x;
                    Rsu_pos_y(count2) = Rsu.pos_y;
                    count2 = count2 + 1; 
                end
            end
            h_rsu_o = plot(Rsu_pos_x, Rsu_pos_y, "dr");

            xlim([0 1000]); ylim([0 1000]);
            xlabel("X Position (meters)");
            ylabel("Y Position (meters)")
            handles = {h_gNB_Bs_o, h_eMBB_UE_o, h_mMTC_UE_o, h_veh_UE_cellular_l2r_o, h_veh_UE_cellular_r2l_o, h_veh_UE_sidelink_l2r_o, h_veh_UE_sidelink_r2l_o, h_rsu_o};
            names = ["BS", "eMBB-UE", "mMTC-UE", "Vehicular-UE-Cellular(l2r)", "Vehicular-UE-Celluar(r2l)", "Vehicular-UE-Sidelink(l2r)", "Vehicular-UE-Sidelink(r2l)", "Rsu"];
            names = names(~cellfun('isempty', handles));
            handles = handles(~cellfun('isempty', handles));
            legend([handles{:}], names, location="northeastoutside");
        end            

        function [] = run_winner2_model_simulation_drop(self, t_i)
            % run the WINNER2 channel model in every simulation drop, in
            % order to get the channel gains and the SINRs in every links.
            % The following codes are similar to the following example,
            % Example: https://ww2.mathworks.cn/help/comm/ug/simultaneous-simulation-of-multiple-fading-channels-with-winner-ii-channel-model.html
            % The usage of the winner2 model can refer to
            % https://ww2.mathworks.cn/help/comm/ref/winner2.layoutparset.html
            % https://ww2.mathworks.cn/help/comm/ref/comm.winner2channel-system-object.html
            % https://ww2.mathworks.cn/help/comm/ref/winner2.wimparset.html
            % https://ww2.mathworks.cn/help/comm/ug/mapping-winner-ii-public-download-to-winner2channel.html		
            AA(1) = self.ULA4_Bs;
            AA(2) = self.ULA2_Ue_Car;
            AA_INDEX_ONE = 1;
            AA_INDEX_TWO = 2;
            
            % decide which vehicular UEs are using UL cellular, and which
            % vehicular UEs are using sidelink
            self.veh_UEs_in_cellular = cell(self.M_veh, 1);
            self.veh_UEs_in_sidelink = cell(self.M_veh, 1);
            next_UE_in_cellular = 1;
            next_UE_in_sidelink = 1;
            for count = 1 : self.M_veh
                if self.pd_veh_celluar_or_sidelink.random() == 0
                    self.veh_UEs_in_cellular{next_UE_in_cellular} = self.vehicular_UEs{count};
                    next_UE_in_cellular = next_UE_in_cellular + 1;
                else
                    self.veh_UEs_in_sidelink{next_UE_in_sidelink} = self.vehicular_UEs{count};
                    next_UE_in_sidelink = next_UE_in_sidelink + 1;
                end
            end
            self.veh_UEs_in_cellular = self.veh_UEs_in_cellular(~cellfun('isempty', self.veh_UEs_in_cellular));
            self.veh_UEs_in_sidelink = self.veh_UEs_in_sidelink(~cellfun('isempty', self.veh_UEs_in_sidelink));
            self.M_veh_UE_cellular = numel(self.veh_UEs_in_cellular);
            self.M_veh_UE_sidelink = numel(self.veh_UEs_in_sidelink);

            % network characteristics for the vehicular UEs either in
            % cellular mode or in sidelink mode
            % The effective size will be dynamically determined by
            % M_veh_UE_cellular, M_veh_UE_sidelink in every simulation drop
            self.nw_characteristics_v2x_CL.SINRdB{t_i} = cell(self.M_veh_UE_cellular, 1);
            self.nw_characteristics_v2x_CL.nr_packets{t_i} = cell(self.M_veh_UE_cellular, 1);

            self.nw_characteristics_v2x_SL.SINRdB{t_i} = cell(self.M_veh_UE_sidelink, 1);
            self.nw_characteristics_v2x_SL.nr_packets{t_i} = cell(self.M_veh_UE_sidelink, 1);               

            % refer to
            % https://ww2.mathworks.cn/help/comm/ref/winner2.layoutparset.html
            % for the layoutparseet configuration
            TxAntArrayIdx = num2cell([AA_INDEX_ONE; ...
                AA_INDEX_TWO * ones(self.M_veh_UE_sidelink, 1)]);   % Index in antenna array inventory vector of gNB and sidelinke 
            RxAntArrayIdx = AA_INDEX_TWO * ones(1, self.M_eMBB + self.M_mMTC + self.M_veh_UE_cellular + self.M_rsu); %  Index in antenna array inventory vector for the single gNB and all the RSUs
            numLinks = self.M_eMBB + self.M_mMTC + self.M_veh_UE_cellular + self.nr_corners_per_cluster * self.M_veh_UE_sidelink; 
            range = Constants_c.diameter; % Layout range (meters)
            cfgLayout = winner2.layoutparset(RxAntArrayIdx, TxAntArrayIdx, numLinks, AA, range);
 
            % refer to the report of the description of the arrange of the
            % links in the ul network model 
            firstTx_Bs = 1;
            firstTx_veh_UE_sidelink = firstTx_Bs + 1;
            lastTx_veh_UE_sidelink = firstTx_Bs + self.M_veh_UE_sidelink;

            first_eMBB_UE_sindex = lastTx_veh_UE_sidelink + 1;
            last_eMBB_UE_sindex = lastTx_veh_UE_sidelink + self.M_eMBB;
            
            first_mMTC_UE_sindex = last_eMBB_UE_sindex + 1;
            last_mMTC_UE_sindex = last_eMBB_UE_sindex + self.M_mMTC;

            first_veh_UE_cellular_sindex = last_mMTC_UE_sindex + 1;
            last_veh_UE_cellular_sindex = last_mMTC_UE_sindex + self.M_veh_UE_cellular;
            
            pairing_Bs_eMBB_Ue = [ones(1, self.M_eMBB); first_eMBB_UE_sindex : last_eMBB_UE_sindex];
            pairing_Bs_mMTC_Ue = [ones(1, self.M_mMTC); first_mMTC_UE_sindex : last_mMTC_UE_sindex];
            pairing_Bs_veh_cellular_Ue = [ones(1, self.M_veh_UE_cellular); first_veh_UE_cellular_sindex : last_veh_UE_cellular_sindex];

            if (4 ~= Constants_c.nr_clusters)
                error('only number of 4 clusters are supported.')
            end
            rsu_cluster_first_sindex = zeros(Constants_c.nr_clusters, 1);
            rsu_cluster_last_sindex = zeros(Constants_c.nr_clusters, 1);
            for count = 1 : Constants_c.nr_clusters
                rsu_cluster_first_sindex(count) = last_veh_UE_cellular_sindex + (count - 1) * self.nr_corners_per_cluster + 1;
                rsu_cluster_last_sindex(count) = last_veh_UE_cellular_sindex + (count) * self.nr_corners_per_cluster;
            end
            
            pairing_sidelink_UE_rsu = zeros(2, self.nr_corners_per_cluster * self.M_veh_UE_sidelink);           
            for count = 1 : self.M_veh_UE_sidelink
                veh_ue_sindex = firstTx_veh_UE_sidelink + count - 1;
                cluster = self.veh_UEs_in_sidelink{count}.cluster;
                first_rsu_sindex = rsu_cluster_first_sindex(cluster);
                last_rsu_sindex = rsu_cluster_last_sindex(cluster);
                pairing_cluster = [veh_ue_sindex * ones(1, self.nr_corners_per_cluster); first_rsu_sindex: last_rsu_sindex];
                pairing_sidelink_UE_rsu(:, ...
                    (1 + (count - 1) * self.nr_corners_per_cluster) : (count * self.nr_corners_per_cluster)) = pairing_cluster;
            end

            
            cfgLayout.Pairing = [pairing_Bs_eMBB_Ue pairing_Bs_mMTC_Ue pairing_Bs_veh_cellular_Ue pairing_sidelink_UE_rsu];   
            WINNER2_B1 = 3;  % 3 for B1 in winner2 model
            cfgLayout.ScenarioVector = WINNER2_B1 * ones(1, numLinks); 
%             LOS = 1; % 0 for NLOS; 1 for LOS
%             cfgLayout.PropagConditionVector = LOS * ones(1, numLinks); 
            % setup the positions
            cfgLayout.Stations(1).Pos(1:2) = [self.gNB_Bs.pos_x; self.gNB_Bs.pos_y];
            % setup the sidelink TX UE
            for sindex = firstTx_veh_UE_sidelink : lastTx_veh_UE_sidelink
                x = round(self.veh_UEs_in_sidelink{sindex - firstTx_veh_UE_sidelink + 1}.pos_x);
                y = round(self.veh_UEs_in_sidelink{sindex - firstTx_veh_UE_sidelink + 1}.pos_y);
                cfgLayout.Stations(sindex).Pos(1:2) = [x; y];      
                cfgLayout.Stations(sindex).Velocity = zeros(3, 1);
                if (strcmpi(self.veh_UEs_in_sidelink{sindex - firstTx_veh_UE_sidelink + 1}.direction, 'l2r'))
                    cfgLayout.Stations(sindex).Velocity(1) = Constants_c.velocity;
                elseif (strcmpi(self.veh_UEs_in_sidelink{sindex - firstTx_veh_UE_sidelink + 1}.direction, 'r2l'))
                    cfgLayout.Stations(sindex).Velocity(1) = -Constants_c.velocity;
                else
                    error('either l2r or r2l is allowed!')
                end                
            end

            for sindex = first_eMBB_UE_sindex : last_eMBB_UE_sindex
                x = round(self.eMBB_UEs{sindex - first_eMBB_UE_sindex + 1}.pos_x);
                y = round(self.eMBB_UEs{sindex - first_eMBB_UE_sindex + 1}.pos_y);
                cfgLayout.Stations(sindex).Pos(1:2) = [x; y];      
                cfgLayout.Stations(sindex).Velocity = rand(3,1) - 0.5;
            end
            for sindex = first_mMTC_UE_sindex : last_mMTC_UE_sindex
                x = round(self.mMTC_UEs{sindex - first_mMTC_UE_sindex + 1}.pos_x);
                y = round(self.mMTC_UEs{sindex - first_mMTC_UE_sindex + 1}.pos_y);
                cfgLayout.Stations(sindex).Pos(1:2) = [x; y];    
                cfgLayout.Stations(sindex).Velocity = 1e-10 * (rand(3,1) - 0.5); % very very small veloctiy
            end
            for sindex = first_veh_UE_cellular_sindex : last_veh_UE_cellular_sindex
                x = round(self.veh_UEs_in_cellular{sindex - first_veh_UE_cellular_sindex + 1}.pos_x);
                y = round(self.veh_UEs_in_cellular{sindex - first_veh_UE_cellular_sindex + 1}.pos_y);
                cfgLayout.Stations(sindex).Pos(1:2) = [x; y];    
                cfgLayout.Stations(sindex).Velocity = zeros(3,1);
                if (strcmpi(self.veh_UEs_in_cellular{sindex - first_veh_UE_cellular_sindex + 1}.direction, 'l2r'))
                    cfgLayout.Stations(sindex).Velocity(1) = Constants_c.velocity;
                elseif (strcmpi(self.veh_UEs_in_cellular{sindex - first_veh_UE_cellular_sindex + 1}.direction, 'r2l'))
                    cfgLayout.Stations(sindex).Velocity(1) = -Constants_c.velocity;
                else
                    error('either l2r or r2l is allowed!')
                end
            end
            for count = 1 : Constants_c.nr_clusters
                start_sindex = rsu_cluster_first_sindex(count);
                last_sindex = rsu_cluster_last_sindex(count);
                for sindex = start_sindex : last_sindex
                    x = round(self.Rsu_clusters{count}(sindex - start_sindex + 1).pos_x);
                    y = round(self.Rsu_clusters{count}(sindex - start_sindex + 1).pos_y);
                    cfgLayout.Stations(sindex).Pos(1:2) = [x; y]; 
                    cfgLayout.Stations(sindex).Velocity = 1e-10 * (rand(3,1) - 0.5); % very very small veloctiy
                end
            end

			% refer to
            % https://ww2.mathworks.cn/help/comm/ref/winner2.wimparset.html
            % for the wimparset configuration
            cfgWim = winner2.wimparset;
            cfgWim.NumTimeSamples = self.winner2_channel_smp_time;
            cfgWim.IntraClusterDsUsed = "yes";
            cfgWim.CenterFrequency = 2.6e9;
            cfgWim.UniformTimeSampling = "no";
            cfgWim.ShadowingModelUsed = "yes";
            cfgWim.PathLossModelUsed = "yes";
            cfgWim.UseManualPropCondition = "no";
            cfgWim.RandomSeed = 31415926 + t_i; % For repeatability
			
            % refer to
            % https://ww2.mathworks.cn/help/comm/ref/comm.winner2channel-system-object.html
            % for the winner2channel configuration           
            self.WINNERChan = comm.WINNER2Channel(cfgWim, cfgLayout);
            chanInfo = info(self.WINNERChan);
            frameLen = self.winner2_channel_smp_time;
            txSig = cellfun(@(x) [ones(1, x); zeros(frameLen - 1, x)], ...
            num2cell(chanInfo.NumBSElements)', UniformOutput = false);
            [~, path_gains] = self.WINNERChan(txSig); % Pass impulse signal through each link
			
			% calculate the SINR for every links
            next_link = 1;
            Pcmax_dBm = 23; % Assuming UE's maximum transmission power is 23dBm
            P_RB_dBm = Pcmax_dBm - 10 * log10(Constants_c.N_RB);
            P_noise_RB_dBm = -174 + 10 * log10(Constants_c.f_sc * Constants_c.N_sc_RB);
            rx_antenna_gain = 5; % antenna gain in UL is 5dB in the gNB (in Table1 in paper)
            for count = 1 : self.M_eMBB
                path_gains_of_link = path_gains{next_link};
                next_link = next_link + 1;
                mean_path_gain_dB = 10 * log10(mean(abs(path_gains_of_link) .^ 2, 'all'));
                SINR_dB = P_RB_dBm + mean_path_gain_dB - P_noise_RB_dBm + rx_antenna_gain;
                self.eMBB_UEs{count}.sinr_dB = SINR_dB;
            end
            for count = 1 : self.M_mMTC
                path_gains_of_link = path_gains{next_link};
                next_link = next_link + 1;
                mean_path_gain_dB = 10 * log10(mean(abs(path_gains_of_link) .^ 2, 'all'));
                SINR_dB = P_RB_dBm + mean_path_gain_dB - P_noise_RB_dBm + rx_antenna_gain;
                self.mMTC_UEs{count}.sinr_dB = SINR_dB;
            end 
            for count = 1 : self.M_veh_UE_cellular
                path_gains_of_link = path_gains{next_link};
                next_link = next_link + 1;
                mean_path_gain_dB = 10 * log10(mean(abs(path_gains_of_link) .^ 2, 'all'));
                SINR_dB = P_RB_dBm + mean_path_gain_dB - P_noise_RB_dBm + rx_antenna_gain;
                self.veh_UEs_in_cellular{count}.sinr_dB = SINR_dB;
            end 
			% Special handling in the calculation of the SINR for sidelink UE, where
			% the RSUs in the corners of the cluster will be used to measure the SINR
            start_sidelink_link = next_link;
            for count = 1 : self.M_veh_UE_sidelink
                rsu_links = (start_sidelink_link: start_sidelink_link + self.nr_corners_per_cluster - 1) + (count - 1) * self.nr_corners_per_cluster;
                sinr_dB_rsu = zeros(self.nr_corners_per_cluster, 1);
                for link = rsu_links
                    path_gains_of_link = path_gains{link};
                    mean_path_gain_dB = 10 * log10(mean(abs(path_gains_of_link) .^ 2, 'all'));
                    SINR_dB = P_RB_dBm + mean_path_gain_dB - P_noise_RB_dBm;
                    sinr_dB_rsu(link - rsu_links(1) + 1) = SINR_dB;
                end
                self.veh_UEs_in_sidelink{count}.sinr_dB = mean(sinr_dB_rsu);
            end
        end

        function [prb_utilization_v2x_o, ...
                prb_utilization_eMBB_o, ...
                prb_utilization_mMTC_o, ...
                prb_utilization_o, ...
                outage_rate_v2x_o, ...
                outage_rate_eMBB_o, ...
                outage_rate_mMTC_o, ...
                outage_rate_o] = start_simulation(self)
            self.init_simulation();
            self.init_winnerII_model();
            self.init_slicing_ratio(1/3, 1/3, 1/3);
            for r = 1 : self.sim_nr_round
                self.prepare_simulation_round(r)
                for t = 1 : self.sim_nr_drops_per_round
                    fprintf(1, '.');
                    self.prepare_simulation_drop(t);                
                    self.run_winner2_model_simulation_drop(t);
                    self.collect_nw_characteristic(t);
                    self.calculate_PRB_utilization(t, 'ul');            
                end
                [self.avg_prb_utilization_ratio_v2x(r), ...
                    self.avg_prb_utilization_ratio_eMBB(r), ...
                    self.avg_prb_utilization_ratio_mMTC(r), ...
                    self.avg_prb_utilization_ratio(r)] = self.calculate_avg_PRB_utilization();  
                [self.avg_outage_rate_v2x(r), ...
                    self.avg_outage_rate_eMBB(r), ...
                    self.avg_outage_rate_mMTC(r), ...
                    self.avg_outage_rate(r)] = self.calculate_avg_outage_rate();                
                fprintf(Constants_c.fid, '\n');
                fprintf(Constants_c.fid, 'avg prb utilization at round (%d) = %.3f\n', r, self.avg_prb_utilization_ratio(r));
                % skip the update of the slicing ratio in the last round
                if r == self.sim_nr_round
                    break;
                else
                    self.update_slicing_ratio(r);
                end
            end   
            fprintf(Constants_c.fid, 'Summary:\n');
            for r = 1 : self.sim_nr_round
                fprintf(Constants_c.fid, 'avg prb utilization at round (%d) %.3f %.3f %.3f %.3f\n', r, ...
                    self.avg_prb_utilization_ratio_v2x(r), ...
                    self.avg_prb_utilization_ratio_eMBB(r), ...
                    self.avg_prb_utilization_ratio_mMTC(r), ...
                    self.avg_prb_utilization_ratio(r));
            end

            for r = 1 : self.sim_nr_round
                fprintf(Constants_c.fid, 'avg outage at round (%d) %.3f %.3f %.3f %.3f\n', r, ...
                    self.avg_outage_rate_v2x(r), ...
                    self.avg_outage_rate_eMBB(r), ...
                    self.avg_outage_rate_mMTC(r), ...
                    self.avg_outage_rate(r));
            end

            prb_utilization_v2x_o = mean(self.avg_prb_utilization_ratio_v2x);
            prb_utilization_eMBB_o = mean(self.avg_prb_utilization_ratio_eMBB);
            prb_utilization_mMTC_o = mean(self.avg_prb_utilization_ratio_mMTC);
            prb_utilization_o = mean(self.avg_prb_utilization_ratio);

            outage_rate_v2x_o = mean(self.avg_outage_rate_v2x);
            outage_rate_eMBB_o = mean(self.avg_outage_rate_eMBB);
            outage_rate_mMTC_o = mean(self.avg_outage_rate_mMTC);
            outage_rate_o = mean(self.avg_outage_rate);
        end         
    end

end