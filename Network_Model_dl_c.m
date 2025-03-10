classdef Network_Model_dl_c < Network_Model_c
    properties
    end

    methods
        function self = Network_Model_dl_c(M_eMBB_i, M_veh_i, M_mMTC_i, lambda_e_i, sim_time_i)
            self@Network_Model_c(M_eMBB_i, M_veh_i, M_mMTC_i, lambda_e_i, sim_time_i);
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
            numTx = 1; % single gNB
            numRx = length(RxAntArrayIdx);   
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
            [rxSig, path_gains] = self.WINNERChan(txSig); % Pass impulse signal through each link
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
    end
end