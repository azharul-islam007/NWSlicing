classdef RL_heuristic_dl_c < slicing_ratio_calculator_c
    properties
        pd_v2x_SINR_dB_;
        pd_v2x_Nr_packet_;

        pd_eMBB_SINR_dB_; 
        pd_eMBB_nr_session_;

        pd_mMTC_SINR_dB_;
        pd_mMTC_log2_packet_size_;
        pd_mMTC_tx_;
    end

    methods
        function self = RL_heuristic_dl_c(alpha_i, pd_v2x_SINR_dB_i, pd_v2x_Nr_packet_i, ...
                pd_eMBB_SINR_dB_i, pd_eMBB_nr_session_i, ...
                pd_mMTC_SINR_dB_i, pd_mMTC_log2_packet_size_i, pd_mMTC_tx_i)
             self@slicing_ratio_calculator_c(alpha_i);

             self.pd_v2x_SINR_dB_ = pd_v2x_SINR_dB_i;
             self.pd_v2x_Nr_packet_ = pd_v2x_Nr_packet_i;

             self.pd_eMBB_SINR_dB_ = pd_eMBB_SINR_dB_i;
             self.pd_eMBB_nr_session_ = pd_eMBB_nr_session_i;

             self.pd_mMTC_SINR_dB_ = pd_mMTC_SINR_dB_i;
             self.pd_mMTC_log2_packet_size_ = pd_mMTC_log2_packet_size_i;
             self.pd_mMTC_tx_ = pd_mMTC_tx_i;
        end

        function [alpha_sel_o] = offline_RL_heuristic_learning_slicing_ratio(self, rand_seed_i, ...
                nr_episode_i, tao_i, Avg_T_i, V_v2x_i, M_eMBB_i, M_mMTC_i)
            rand_seed = rand_seed_i;
            % The temperature parameter in Table 1 in the paper, refer to
            % section 'SELECTION CRITERION' for its impact on the selection
            % probability
            tao = tao_i;     
            % number of episodes. Every 'action' is evaluated in one
            % episode and a new 'action' will be taken after the
            % completition of an episode.
            nr_episode = nr_episode_i;  
            % Averaging window before taking the next action
            Avg_T = Avg_T_i;  
            % The number of vehicular UE in each cluster
            V_v2x = V_v2x_i;
            % The number of eMBB UE
            M_eMBB = M_eMBB_i;
            % The number of mMTC UE
            M_mMTC = M_mMTC_i;

            rng(rand_seed);   % seed initialization for repeatibility
            
            B = Constants_c.B;     
            Fd = Constants_c.Fd;
            Sm = Constants_c.Sm;  
            N_DL = Constants_c.N_RB; 
            Nr_S = Constants_c.Nr_S;                   
            Ar = Constants_c.Ar;                     
            Ax = Constants_c.Ax;                  
            V2X = Constants_c.I_V2X;
            eMBB = Constants_c.I_eMBB;
            mMTC = Constants_c.I_mMTC;

            alpha = self.alpha_;              
            beta_DL = zeros(Nr_S, Ar, Ax);  % Actions of the RL algorithm

            for r = 1 : Ar
                for k = 1 : Ax
                    % slicing ratio between V2X and eMBB service, which is the action
                    % as the result from the RL algorithm for slice 1 and slice 2
                    beta_DL(V2X, r, k) = (1 - 0.05 * r) * (0.05 * k);
                    beta_DL(eMBB, r, k) = (1 - 0.05 * r) * (1 - 0.05 * k);
                    % slicing ratio for mMTC in slice 3
                    beta_DL(mMTC, r, k) = 0.05 * r;
                end
            end
            
            R_DL = cell(Nr_S, 1);       % reward after taking an action
            C = Constants_c.nr_clusters; % number of cluster

            % The off-line network model for the V2X serive, which should simulate the
            % characteristics of the real network model
            pd_v2x_SINR_dB = cell(C, 1);
            pd_v2x_Nr_packet = cell(C, 1);
            for j = 1 : C
                Vj = V_v2x(j);
                if (Vj > 0)
                    pd_v2x_SINR_dB{j} = cell(Vj, 1);
                    pd_v2x_Nr_packet{j} = cell(Vj, 1);
                    for i = 1 : Vj
                        pd_v2x_SINR_dB{j}{i} = self.pd_v2x_SINR_dB_; 
                        pd_v2x_Nr_packet{j}{i} = self.pd_v2x_Nr_packet_; 
                    end
                end
            end

            % The off-line network model for the eMBB serive, which should simulate the
            % characteristics of the real network model         
            pd_eMBB_nr_session = cell(M_eMBB, 1);
            pd_eMBB_SINR_dB = cell(M_eMBB, 1);
            for m = 1 : M_eMBB
                pd_eMBB_SINR_dB{m} = self.pd_eMBB_SINR_dB_; 
                pd_eMBB_nr_session{m} = self.pd_eMBB_nr_session_; 
            end

            % The off-line network model for the mMTC service, which should simulate
            % the characteristics of the real network model          
            pd_mMTC_log2_packet_size = cell(M_mMTC, 1);
            pd_mMTC_tx = cell(M_mMTC, 1);
            pd_mMTC_SINR_dB = cell(M_mMTC, 1);
            for m = 1 : M_mMTC
                pd_mMTC_log2_packet_size{m} = self.pd_mMTC_log2_packet_size_; 
                pd_mMTC_tx{m} = self.pd_mMTC_tx_; 
                pd_mMTC_SINR_dB{m} = self.pd_mMTC_SINR_dB_;
            end

            Q_DL = zeros(Ar, Ax, nr_episode);        % Q-table
            actions_r = zeros(nr_episode, 1);
            actions_x = zeros(nr_episode, 1);
            Psi_DL = cell(Nr_S, nr_episode); % normalized resource utilization after taking an action, for three services

            for episode = 1 : nr_episode
                % [1] apply the action into the off-line network model
                % [1.1] For V2X service, to get the number of arriving packets and the
                % SINR for the vehicular UEs during the averaging window period from
                % the off-line network model
                [m_v2x_DL, Speff_v2x_DL] = offline_nw_v2x_DL_f(...
                    Avg_T, C, V_v2x, pd_v2x_SINR_dB, pd_v2x_Nr_packet);

                % [1.2] For eMBB serice, to get the number of arriving sessions and the
                % SINR for the vehicular UEs during the averaging window period from
                % the off-line network model
                [Rb, Speff_eMBB] = offline_nw_eMBB_DL_f(Avg_T, ...
                    M_eMBB, ...
                    pd_eMBB_nr_session, ...
                    pd_eMBB_SINR_dB);

                % [1.3] For mMTC serice, to get the number of sporadic transmissions,
                % the packet size and the SINR for the mMTC UEs during the
                % averaging window period from the off-line network model
                [tx_mMTC_DL, ps_mMTC_DL, Speff_mMTC_DL] = offline_nw_mMTC_DL_f(...
                    Avg_T, ...
                    M_mMTC, ...
                    pd_mMTC_SINR_dB, ...
                    pd_mMTC_log2_packet_size, ...
                    pd_mMTC_tx);

                % [1.3] Action selection based on the Q-table and the softmax policy
                % Refer to the report for the detail
                if (episode == 1)
                    exp_Qr = ones(Ar, 1);
                    exp_Qx = ones(Ax, 1);
                else
                    exp_Q = exp(Q_DL(:, :, episode - 1) / tao);
                    exp_Qr = sum(exp_Q, 2);
                    exp_Qx = sum(exp_Q, 1);
                end
                pr = exp_Qr / sum(exp_Qr);
                px = exp_Qx / sum(exp_Qx);
                pd_r = makedist('Multinomial', 'Probabilities', pr);
                pd_x = makedist('Multinomial', 'Probabilities', px);
                ar = pd_r.random();
                ax = pd_x.random();
                if ar < 1 || ar > Ar || ax < 1 || ax > Ax
                    error('action should be from 1 to Ar, or 1 to Ax');
                end
                actions_r(episode) = ar;
                actions_x(episode) = ax;

                % [1.4] Apply the action to the off-line network and calculate the
                % resource utilization.
                % For V2X service, calculate the resource
                % utilization according to equation [3] in the paper
                m = m_v2x_DL;
                Speff_v2x = Speff_v2x_DL;
                assert(numel(m) == Avg_T, 'The averaging period for V2X before taking the next action is not consistent.');
                nominator = 0;
                denominator = Avg_T * B * Fd;
                for t = 1: Avg_T
                    assert(numel(m{t}) == C, 'The number of clusters is not consistent');
                    for j = 1: C
                        Vj = numel(m{t}{j});
                        for i = 1 : Vj
                            if (Speff_v2x{t}{j}{i} > 0)
                                nominator = nominator + m{t}{j}{i} * Sm / Speff_v2x{t}{j}{i};
                            end
                        end
                    end
                end
                Gamma1 = nominator / denominator;

                % For eMBB service, calculate the resource
                % utilization according to equation [4] in the paper
                assert(numel(Rb) == Avg_T, 'The number of period for eMBB actions is not cosnsistent');
                nominator = 0;
                denominator = Avg_T;
                for t = 1 : Avg_T
                    M = numel(Rb{t});
                    for m = 1 : M
                        if (Speff_eMBB{t}{m} > 0)
                            rho_x_m_t = Rb{t}{m} / (Speff_eMBB{t}{m} * B);
                            nominator = nominator + rho_x_m_t;
                        end
                    end
                end
                Gamma2 = nominator / denominator;

                % For mMTC service, calculate the resource
                % utilization according to equation [x] in the report
                tx = tx_mMTC_DL;
                ps = ps_mMTC_DL;
                Speff_mMTC = Speff_mMTC_DL;
                assert(numel(tx) == Avg_T, 'The averaging period for mMTC before taking the next action is not consistent.');
                nominator = 0;
                denominator = Avg_T * B * Fd;
                for t = 1: Avg_T
                    assert(numel(tx{t}) == M_mMTC, 'The number of mMTC users is not consistent');
                    for m = 1 : M_mMTC
                        if (Speff_mMTC{t}{m} > 0)
                            nominator = nominator + tx{t}{m} * ps{t}{m} / Speff_mMTC{t}{m};
                        end
                    end
                end
                Gamma3 = nominator / denominator;

                % [2] Reward computation
                % [2.1] Calculate the normalized resource utilization according to equation [10], [12]
                beta_V2X_DL = beta_DL(V2X, ar, ax);
                beta_eMBB_DL = beta_DL(eMBB, ar, ax);
                beta_mMTC_DL = beta_DL(mMTC, ar, ax);
                Psi_DL{V2X, episode} = Gamma1 / (beta_V2X_DL * N_DL);
                Psi_DL{eMBB, episode} = Gamma2 / (beta_eMBB_DL * N_DL);
                Psi_DL{mMTC, episode} = Gamma3 / (beta_mMTC_DL * N_DL);
                % [2.2] Calculate the reward according to equation [13][14]
                R_DL{V2X} = 0;
                R_DL{eMBB} = 0;
                R_DL{mMTC} = 0;
                if Psi_DL{V2X, episode} <= 1
                    R_DL{V2X} = exp(Psi_DL{V2X, episode});
                else
                    R_DL{V2X} = 1/exp(Psi_DL{V2X, episode}); % TODO: 1/Psi_DL{V2X} or 1/exp(Psi_UL{V2X})?
                end
                if Psi_DL{eMBB, episode} <= 1
                    R_DL{eMBB} = exp(Psi_DL{eMBB, episode});
                else
                    R_DL{eMBB} = 1/exp(Psi_DL{eMBB, episode}); % TODO: 1/Psi_DL{eMBB} or 1/exp(Psi_UL{eMBB})?
                end
                if Psi_DL{mMTC, episode} <= 1
                    R_DL{mMTC} = exp(Psi_DL{mMTC, episode});
                else
                    R_DL{mMTC} = 1/exp(Psi_DL{mMTC, episode}); % TODO: 1/Psi_DL{mMTC} or 1/exp(Psi_UL{mMTC})?
                end
                R_TOT_DL = (R_DL{V2X} * R_DL{eMBB} * R_DL{mMTC}) ^ (1/3);

                % [3] Q-learning table update according to equation [15]
                if episode > 1
                    for ri = 1 : Ar
                        for xi = 1 : Ax
                            if (ri == ar && xi == ax)
                                Q_DL(ar, ax, episode) = (1 - alpha) * Q_DL(ar, ax, episode-1) + alpha * R_TOT_DL;
                            else
                                Q_DL(ri, xi, episode) = Q_DL(ri, xi, episode-1);
                            end
                        end
                    end
                else
                    Q_DL(ar, ax, episode) = alpha * R_TOT_DL;
                end
%                 fprintf(Constants_c.fid, 'episode=%04d: ar=%02d, ax=%02d, Gamma1=%3.2f, Gamma2=%3.2f, Gamma3=%3.2f, R_DL(V2X)=%3.2f, R_DL(eMBB)=%3.2f, R_DL(mMTC)=%3.2f, R_TOT_DL=%3.2f, beta_V2X=%3.2f, beta_eMBB=%3.2f, beta_mMTC=%3.2f\n', ...
%                     episode, ar, ax, Gamma1, Gamma2, Gamma3, R_DL{V2X}, R_DL{eMBB}, R_DL{mMTC}, R_TOT_DL, beta_V2X_DL, beta_eMBB_DL, beta_mMTC_DL);
            end

%             % To plot the changes of the Q-table during the RL algorithm.
%             % For simplicity, only 10 instance of the Q-table during the RL are plotted.
%             nr_line_to_plot = 1;
%             episode_to_check = nr_episode : -nr_episode / nr_line_to_plot : 1;
%             figure(1);
%             hold on;
%             for episode = (episode_to_check)
%                 mesh(Q_DL(:, :, episode), 'DisplayName', sprintf('episode=%d', episode));
%             end
%             xlabel('actions ax(eMBB, V2X')
%             ylabel('actions ar(mMTC)');
%             zlabel('Q-values');
%             legend;
%             grid on;
%             hold off;

            % Select the final action (r, k) by selecting the one with the maximum Q
            % value
            Q_final = Q_DL(:, :, nr_episode);
            [~, I] = max(Q_final, [], 'all');
            ar_sel = mod((I - 1), Ar) + 1;
            ax_sel = ceil(I / 4);
            
            % Low complexity heuristic tuning 
            alpha_sel = heuristic_tuning_f(ar_sel, ax_sel, actions_r, actions_x, Psi_DL, beta_DL(:, ar, ax));
            alpha_sel_o = alpha_sel;
        end
    end
end