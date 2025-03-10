classdef Network_Model_dl_c < Network_Model_c
    properties
        sim_nr_round;
        avg_prb_utilization_ratio;
    end

    methods
        function self = Network_Model_dl_c(M_eMBB_i, M_veh_i, M_mMTC_i, lambda_e_i, sim_nr_drops_per_round_i, sim_nr_round_i)
            self@Network_Model_c(M_eMBB_i, M_veh_i, M_mMTC_i, lambda_e_i, sim_nr_drops_per_round_i);
            self.sim_nr_round = sim_nr_round_i;
            self.avg_prb_utilization_ratio = zeros(self.sim_nr_round, 1);
        end

        function [] = update_slicing_ratio(self, round_i)
            if (Constants_c.ignore_online_nw_characteristics)
                pd_v2x_SINR_dB = Constants_c.pd_SINR_dB_by_default;
                pd_v2x_Nr_packet = Constants_c.pd_Nr_packet_by_default;
                pd_eMBB_SINR_dB = Constants_c.pd_SINR_dB_by_default;
                pd_eMBB_nr_session = Constants_c.pd_Nr_session_by_default;
                pd_mMTC_SINR_dB = Constants_c.pd_SINR_dB_by_default;
                pd_mMTC_tx = Constants_c.pd_tx_by_default;
                pd_mMTC_log2_packet_size = Constants_c.pd_log2_packet_size_by_default;
            else
                [pd_v2x_SINR_dB, pd_v2x_Nr_packet, ...
                    pd_eMBB_SINR_dB, pd_eMBB_nr_session, ...
                    pd_mMTC_SINR_dB, pd_mMTC_log2_packet_size, pd_mMTC_tx] = self.estimate_nw_parameters_f();
            end

            alpha = Constants_c.alpha_learning_rate; % RL learning rate
            RL_heuristic_dl = RL_heuristic_dl_c(alpha, pd_v2x_SINR_dB, pd_v2x_Nr_packet, ...
                pd_eMBB_SINR_dB, pd_eMBB_nr_session, ...
                pd_mMTC_SINR_dB, pd_mMTC_log2_packet_size, pd_mMTC_tx);

            rand_seed = round_i;
            V_v2x = ones(1, Constants_c.nr_clusters) * floor(self.M_veh / Constants_c.nr_clusters); % Assuming that all Veh UE are uniformally distributed in each cluster
            rem = mod(self.M_veh, Constants_c.nr_clusters);
            % distribute the remaining UEs randomly in any one of the
            % clusters
            if (rem > 0)
                for i = 1 : rem
                    cluster_index = randsample(Constants_c.nr_clusters, 1);
                    V_v2x(cluster_index) = V_v2x(cluster_index) + 1;
                end
            end

            M_eMBB = self.M_eMBB;
            M_mMTC = self.M_mMTC;
            nr_episode = Constants_c.nr_episode;
            tao = Constants_c.tao;
            Avg_T = Constants_c.Avg_T;
            alpha_sel = RL_heuristic_dl.offline_RL_heuristic_learning_slicing_ratio(rand_seed, ...
                nr_episode, tao, Avg_T, V_v2x, M_eMBB, M_mMTC);            
            V2X = Constants_c.I_V2X;
            eMBB = Constants_c.I_eMBB;
            mMTC = Constants_c.I_mMTC;
            self.alpha_v2x = alpha_sel(V2X);
            self.alpha_eMBB = alpha_sel(eMBB);
            self.alpha_mMTC = alpha_sel(mMTC);
            fprintf(1, 'slicing ratio updated after round (%d) = %.3f %.3f %.3f\n', round_i, self.alpha_v2x, self.alpha_eMBB, self.alpha_mMTC);
        end

        function [pd_v2x_SINR_dB_o, pd_v2x_Nr_packet_o, ...
                pd_eMBB_SINR_dB_o, pd_eMBB_nr_session_o, ...
                pd_mMTC_SINR_dB_o, pd_mMTC_log2_packet_size_o, pd_mMTC_tx_o] = estimate_nw_parameters_f(self)
            sinr_dB = reshape(self.nw_characteristics_v2x.SINRdB, [], 1);
            pd_v2x_SINR_dB_o = fitdist(sinr_dB, 'Normal'); 
            nr_packets = reshape(self.nw_characteristics_v2x.nr_packets, [], 1);            
            pd_v2x_Nr_packet_o = fitdist(nr_packets, 'Poisson');

            sinr_dB = reshape(self.nw_characteristics_eMBB.SINRdB, [], 1);
            pd_eMBB_SINR_dB_o = fitdist(sinr_dB, 'Normal'); 
            nr_session = reshape(self.nw_characteristics_eMBB.nr_session, [], 1);            
            pd_eMBB_nr_session_o = fitdist(nr_session, 'Poisson');

            sinr_dB = reshape(self.nw_characteristics_mMTC.SINRdB, [], 1);
            pd_mMTC_SINR_dB_o = fitdist(sinr_dB, 'Normal'); 

            packet_size = reshape(self.nw_characteristics_mMTC.packet_size, [], 1);
            [GC, GR] = groupcounts(packet_size);
            log2_GR = log2(GR);
            max_log2_GR = max(log2_GR);
            prob = zeros(max_log2_GR, 1);
            for count = 1 : length(log2_GR)
                k = log2_GR(count);
                prob(k) = GC(count) / sum(GC);
            end
            pd_mMTC_log2_packet_size_o = makedist('Multinomial', 'Probabilities', prob);

            tx = reshape(self.nw_characteristics_mMTC.tx_flags, [], 1);
            pd_mMTC_tx_o = fitdist(tx, 'Binomial', 'N', 1);           
        end

        function [] = plot(self)
            [h_gNB_Bs_o, h_eMBB_UE_o, h_mMTC_UE_o] = self.plot@Network_Model_c();         
            % special plotting for the vehicular UEs in the DL network
            % model is added here
            nr_veh_UE_l2r = 0;
            nr_veh_UE_r2l = 0;
            for count = 1 : self.M_veh
                if strcmpi(self.vehicular_UEs{count}.direction, 'l2r')
                    nr_veh_UE_l2r = nr_veh_UE_l2r + 1;
                elseif strcmpi(self.vehicular_UEs{count}.direction, 'r2l')
                    nr_veh_UE_r2l = nr_veh_UE_r2l + 1;
                end
            end
            if (nr_veh_UE_l2r > 0)
                veh_UE_pos_x = zeros(nr_veh_UE_l2r, 1);
                veh_UE_pos_y = zeros(nr_veh_UE_l2r, 1);
                count2 = 1;
                for count = 1 : self.M_veh
                    if strcmpi(self.vehicular_UEs{count}.direction, 'l2r')
                        veh_UE_pos_x(count2) = self.vehicular_UEs{count}.pos_x;
                        veh_UE_pos_y(count2) = self.vehicular_UEs{count}.pos_y;
                        count2 = count2 + 1;
                    end
                end
                h_veh_UE_l2r_o = plot(veh_UE_pos_x, veh_UE_pos_y, ">k");
            end
            if (nr_veh_UE_r2l > 0)
                veh_UE_pos_x = zeros(nr_veh_UE_r2l, 1);
                veh_UE_pos_y = zeros(nr_veh_UE_r2l, 1);
                count2 = 1;
                for count = 1 : self.M_veh
                    if strcmpi(self.vehicular_UEs{count}.direction, 'r2l')
                        veh_UE_pos_x(count2) = self.vehicular_UEs{count}.pos_x;
                        veh_UE_pos_y(count2) = self.vehicular_UEs{count}.pos_y;
                        count2 = count2 + 1;
                    end
                end
                h_veh_UE_r2l_o = plot(veh_UE_pos_x, veh_UE_pos_y, "<m");
            end
          
            for count = 1 : self.M_veh
                x = [self.gNB_Bs.pos_x self.vehicular_UEs{count}.pos_x];
                y = [self.gNB_Bs.pos_y self.vehicular_UEs{count}.pos_y];
                plot(x, y, "--r");  
                text(x(2),y(2), sprintf('sinr=%.1f dB', self.vehicular_UEs{count}.sinr_dB));
            end
            xlim([0 1000]); ylim([0 1000]);
            xlabel("X Position (meters)");
            ylabel("Y Position (meters)")
            handles = {h_gNB_Bs_o, h_eMBB_UE_o, h_mMTC_UE_o, h_veh_UE_l2r_o, h_veh_UE_r2l_o};
            names = ["BS","eMBB-UE", "mMTC-UE", "Vehicular-UE(l2r)", "Vehicular-UE(r2l)"];
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

            % refer to
            % https://ww2.mathworks.cn/help/comm/ref/winner2.layoutparset.html
            % for the layoutparseet configuration
            TxAntArrayIdx = {AA_INDEX_ONE}; % Index in antenna array inventory vector for the single gNB
            RxAntArrayIdx = AA_INDEX_TWO * ones(1, self.M_eMBB + self.M_veh + self.M_mMTC); % Index in antenna array inventory vector for the UEs
%             numTx = 1; % single gNB
%             numRx = length(RxAntArrayIdx);   
            numLinks = self.M_eMBB + self.M_veh + self.M_mMTC; % 8 Number of links
            range = Constants_c.diameter; % Layout range (meters)
            cfgLayout = winner2.layoutparset(RxAntArrayIdx, TxAntArrayIdx, numLinks, AA, range);
                        
            % refer to the report of the description of the arrange of the
            % links in the dl network model
            Bs_sindex = 1;
            first_eMBB_UE_sindex = Bs_sindex + 1;
            last_eMBB_UE_sindex = Bs_sindex + self.M_eMBB;
            
            first_mMTC_UE_sindex = last_eMBB_UE_sindex + 1;
            last_mMTC_UE_sindex = last_eMBB_UE_sindex + self.M_mMTC;

            first_veh_UE_sindex = last_mMTC_UE_sindex + 1;
            last_veh_UE_sindex = last_mMTC_UE_sindex + self.M_veh;
            
            pairing_Bs_eMBB_Ue = [ones(1, self.M_eMBB); first_eMBB_UE_sindex : last_eMBB_UE_sindex];
            pairing_Bs_mMTC_Ue = [ones(1, self.M_mMTC); first_mMTC_UE_sindex : last_mMTC_UE_sindex];
            pairing_Bs_veh_Ue = [ones(1, self.M_veh); first_veh_UE_sindex : last_veh_UE_sindex];
            
            cfgLayout.Pairing = [pairing_Bs_eMBB_Ue pairing_Bs_mMTC_Ue pairing_Bs_veh_Ue];          
            WINNER2_B1 = 3;  % 3 for B1 in winner2 model
            cfgLayout.ScenarioVector = WINNER2_B1 * ones(1, numLinks); 
            %LOS = 1; % 0 for NLOS; 1 for LOS
            %cfgLayout.PropagConditionVector = LOS * ones(1, numLinks); 
            % setup the positions for tx, rx points in the winner2 model
            % according to the network model
            cfgLayout.Stations(1).Pos(1:2) = [self.gNB_Bs.pos_x; self.gNB_Bs.pos_y];
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
            for sindex = first_veh_UE_sindex : last_veh_UE_sindex
                x = round(self.vehicular_UEs{sindex - first_veh_UE_sindex + 1}.pos_x);
                y = round(self.vehicular_UEs{sindex - first_veh_UE_sindex + 1}.pos_y);
                cfgLayout.Stations(sindex).Pos(1:2) = [x; y];    
                cfgLayout.Stations(sindex).Velocity = zeros(3,1);
                if (strcmpi(self.vehicular_UEs{sindex - first_veh_UE_sindex + 1}.direction, 'l2r'))
                    cfgLayout.Stations(sindex).Velocity(1) = Constants_c.velocity;
                elseif (strcmpi(self.vehicular_UEs{sindex - first_veh_UE_sindex + 1}.direction, 'r2l'))
                    cfgLayout.Stations(sindex).Velocity(1) = -Constants_c.velocity;
                else
                    error('either l2r or r2l is allowed!')
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
            Pcmax = 10e3; % assuming the max transmission power from the gNB is 10w
            Pcmax_dBm = 10*log10(Pcmax);
            P_RB_dBm = Pcmax_dBm - 10 * log10(Constants_c.N_RB); % power is uniformly distributed in every PRB
            P_noise_RB_dBm = -174 + 10 * log10(Constants_c.f_sc * Constants_c.N_sc_RB); % -174dBm/Hz is the white noise power spectrum at T=25 cel. degree.
            for count = 1 : self.M_eMBB
                path_gains_of_link = path_gains{next_link};
                next_link = next_link + 1;
                mean_path_gain_dB = 10 * log10(mean(abs(path_gains_of_link) .^ 2, 'all'));
                SINR_dB = P_RB_dBm + mean_path_gain_dB - P_noise_RB_dBm;
                self.eMBB_UEs{count}.sinr_dB = SINR_dB;
            end
            for count = 1 : self.M_mMTC
                path_gains_of_link = path_gains{next_link};
                next_link = next_link + 1;
                mean_path_gain_dB = 10 * log10(mean(abs(path_gains_of_link) .^ 2, 'all'));
                SINR_dB = P_RB_dBm + mean_path_gain_dB - P_noise_RB_dBm;
                self.mMTC_UEs{count}.sinr_dB = SINR_dB;
            end 
            for count = 1 : self.M_veh
                path_gains_of_link = path_gains{next_link};
                next_link = next_link + 1;
                mean_path_gain_dB = 10 * log10(mean(abs(path_gains_of_link) .^ 2, 'all'));
                SINR_dB = P_RB_dBm + mean_path_gain_dB - P_noise_RB_dBm;
                self.vehicular_UEs{count}.sinr_dB = SINR_dB;
            end              
        end

        function [avg_prb_utilization_o] = start_simulation(self)
            self.init_simulation();
            self.init_winnerII_model();
            self.init_slicing_ratio(1/3, 1/3, 1/3);
            for r = 1 : self.sim_nr_round
                for t = 1 : self.sim_nr_drops_per_round
                    fprintf(1, '.');
                    self.prepare_simulation_drop(t);                
                    self.run_winner2_model_simulation_drop(t);
                    self.collect_nw_characteristic(t);
                    self.calculate_PRB_utilization(t, 'dl');            
                end
                self.avg_prb_utilization_ratio(r) = self.calculate_avg_PRB_utilization();                
                fprintf(1, '\n');
                fprintf(1, 'avg prb utilization at round (%d) %.3f\n', r, self.avg_prb_utilization_ratio(r));
                % skip the update of the slicing ratio in the last round
                if r == self.sim_nr_round
                    break;
                else
                    self.update_slicing_ratio(r);
                end
            end      
            fprintf(1, 'Summary:\n');
            for r = 1 : self.sim_nr_round
                fprintf(1, 'avg prb utilization at round (%d) = %.3f\n', r, self.avg_prb_utilization_ratio(r));
            end
            avg_prb_utilization_o = mean(self.avg_prb_utilization_ratio);
        end        
    end
end