function [alpha_sel_o] = slicing_RL_UL_f(rand_seed_i, nr_episode_i, tao_i, Avg_T_i, V_v2x_UL_i, V_v2x_SL_i, M_eMBB_i, M_mMTC_i)
%% RAN slicing algorithm, UL
% Note that the procedure of RL for UL direction is quite similar to that
% in the UL, so in the following descriptions, only the differences between
% DL and UL is descripted. One can refer to the script, 'slicing_RL_DL_f',
% in the same folder for more descriptions.

rng(rand_seed_i);   
tao = tao_i;  
nr_episode = nr_episode_i; 
B = 30e3 * 12;              
alpha = 0.1;               
Sm = 300;                  
N_UL = 200;                 % Number of PRBs in UL
Nr_S = 3;  
Ar = 4;
Ax = 20;                    
Fd = 0.5e-3;                
beta_UL = zeros(Nr_S, Ar, Ax);  

V2X = 1;
eMBB = 2;
mMTC = 3;
for r = 1 : Ar 
    for k = 1 : Ax
        % slicing ratio between V2X and eMBB service, which is the action
        % as the result from the RL algorithm for slice 1 and slice 2
        beta_UL(V2X, r, k) = (1 - 0.05 * r) * (0.05 * k);   
        beta_UL(eMBB, r, k) = (1 - 0.05 * r) * (1 - 0.05 * k); 
         % slicing ratio for mMTC in slice 3     
        beta_UL(mMTC, r, k) = 0.05 * r;
    end
end

Avg_T = Avg_T_i;  
R_UL = cell(Nr_S, 1);

C = 4;   

% The off-line network model for the V2X serive in UL using the cellular mode 
V_v2x_UL = V_v2x_UL_i;
pd_v2x_UL_SINR_dB = cell(C, 1);
pd_v2x_UL_Nr_packet = cell(C, 1);
for j = 1 : C
    Vj = V_v2x_UL(j);
    pd_v2x_UL_SINR_dB{j} = cell(Vj, 1);
    pd_v2x_UL_Nr_packet{j} = cell(Vj, 1);
    for i = 1 : Vj
        pd_v2x_UL_SINR_dB{j}{i} = makedist('Uniform', 'Lower', -10, 'Upper', 40); 
        pd_v2x_UL_Nr_packet{j}{i} = makedist('Binomial', 'N', 2, 'p', 0.5);
    end
end

% Different from the RL in Downlink, here the off-line model needs to
% simulate the sidelink characteristics (which should also be inputted from the
% on-line network model).
V_v2x_SL = V_v2x_SL_i;
pd_v2x_SL_SINR_dB = cell(C, 1);
pd_v2x_SL_Nr_packet = cell(C, 1);
for j = 1 : C
    Vj = V_v2x_SL(j);
    pd_v2x_SL_SINR_dB{j} = cell(Vj, 1);
    pd_v2x_SL_Nr_packet{j} = cell(Vj, 1);
    for i = 1 : Vj
        pd_v2x_SL_SINR_dB{j}{i} = makedist('Uniform', 'Lower', -10, 'Upper', 40); 
        pd_v2x_SL_Nr_packet{j}{i} = makedist('Binomial', 'N', 2, 'p', 0.5);
    end
end

% The off-line network model for the eMBB serive, which should simulate the
% characteristics of the real network model
M_eMBB = M_eMBB_i;
pd_nr_eMBB_session = cell(M_eMBB, 1);
pd_eMBB_SINR_dB = cell(M_eMBB, 1);
for m = 1 : M_eMBB
    pd_eMBB_SINR_dB{m} = makedist('Uniform', 'Lower', -10, 'Upper', 40); 
    pd_nr_eMBB_session{m} = makedist('Binomial', 'N', 2, 'p', 0.5);
end


% The off-line network model for the mMTC service, which should simulate
% the characteristics of the real network model
% TODO: to integrate the off-line network model with the (simulated)
% on-line network model
M_mMTC = M_mMTC_i;
pd_packet_size_mMTC = cell(M_mMTC, 1);
pd_tx_mMTC = cell(M_mMTC, 1);
pd_mMTC_SINR_dB = cell(M_mMTC, 1);
max_packet_size = 128; % maxumum packet size = 128bytes, the possible packet sizes = 1, 2, 4, 8, 16, 32, 64, 128 bytes
log2_max_packet_size = log2(max_packet_size);
for m = 1 : M_mMTC
    pd_packet_size_mMTC{m} = makedist('Multinomial', 'Probabilities', 1 / log2_max_packet_size * ones(1, log2_max_packet_size));
    pd_tx_mMTC{m} = makedist('Binomial', 'N', 1, 'p', 0.25);
    pd_mMTC_SINR_dB{m} = makedist('Uniform', 'lower', -10, 'upper', 40);
end

Q_UL = zeros(Ar, Ax, nr_episode); 
actions_r = zeros(nr_episode, 1);
actions_x = zeros(nr_episode, 1);
Psi_UL = cell(Nr_S, nr_episode);
for episode = 1 : nr_episode
    % [1] apply the action into the off-line network model, including the
    % sidelink in UL

    % [1.1] For V2X service either in cellular mode or in sidelink mode, to
    % get the number of arriving packets and the SINR for the vehicular UEs
    % during the averaging window period from the off-line network model
    [m_v2x_UL, m_v2x_SL, Speff_v2x_UL, Speff_v2x_SL] = offline_nw_v2x_ULSL_f(...
        Avg_T, C, V_v2x_UL, pd_v2x_UL_SINR_dB, pd_v2x_UL_Nr_packet, ...
        V_v2x_SL, pd_v2x_SL_SINR_dB, pd_v2x_SL_Nr_packet);

    % [1.2] For eMBB serice, to get the number of arriving sessions and the
    % SINR for the vehicular UEs during the averaging window period from
    % the off-line network model    
    [Rb, Speff_eMBB] = offline_nw_eMBB_UL_f(Avg_T, ...
        M_eMBB, ...
        pd_nr_eMBB_session, ...
        pd_eMBB_SINR_dB);

    [tx_mMTC_UL, ps_mMTC_UL, Speff_mMTC_UL] = offline_nw_mMTC_UL_f(...
        Avg_T, ...
        M_mMTC, ...
        pd_mMTC_SINR_dB, ...
        pd_packet_size_mMTC, ...
        pd_tx_mMTC);    

    % [1.3] Action selection based on the Q-table and the softmax policy
    if (episode == 1)
        exp_Qr = ones(Ar, 1);
        exp_Qx = ones(Ax, 1);    
    else
        exp_Q = exp(Q_UL(:, :, episode - 1) / tao);
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

    % For V2X service either in cellular mode or in sidelink mode,
    % calculate the resource utilization according to equation [3] in the
    % paper
    UL = 1;
    SL = 2;
    m = {m_v2x_UL, m_v2x_SL};
    Speff_v2x = {Speff_v2x_UL, Speff_v2x_SL};
    Gamma1 = cell(numel(m), 1);
    for x = [UL, SL]
        assert(numel(m{x}) == Avg_T, 'The averaging period for V2X before taking the next action is not consistent.');
        nominator = 0;
        denominator = Avg_T * B * Fd;      
        for t = 1: Avg_T
            assert(numel(m{x}{t}) == C, 'The number of clusters is not consistent');
            for j = 1: C
                Vj_v2x_x = numel(m{x}{t}{j});
                for i = 1 : Vj_v2x_x
                    if (Speff_v2x{x}{t}{j}{i} > 0)
                        nominator = nominator + m{x}{t}{j}{i} * Sm / Speff_v2x{x}{t}{j}{i};
                    end
                end
            end
        end
        Gamma1{x} = nominator / denominator;
    end  

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
    tx = tx_mMTC_UL;
    ps = ps_mMTC_UL;    
    Speff_mMTC = Speff_mMTC_UL;
    assert(numel(tx) == Avg_T, 'The averaging period for mMTC before taking the next action is not consistent.');
    nominator = 0;
    denominator = Avg_T * B * Fd;      
    for t = 1: Avg_T
        assert(numel(tx{t}) == M_mMTC, 'The number of mMTC users is not consistent');
        for m = 1 : M_mMTC
            nominator = nominator + tx{t}{m} * ps{t}{m} / Speff_mMTC{t}{m};
        end
    end
    Gamma3 = nominator / denominator;       

    % [2] Reward computation

    % [2.1] Calculate the normalized resource utilization according to equation [9], [11]
    beta_V2X_UL = beta_UL(V2X, ar, ax);
    beta_eMBB_UL = beta_UL(eMBB, ar, ax);
    beta_mMTC_UL = beta_UL(mMTC, ar, ax);
    Psi_UL{V2X, episode} = (Gamma1{UL} + Gamma1{SL}) / (beta_V2X_UL * N_UL);
    Psi_UL{eMBB, episode} = Gamma2 / (beta_eMBB_UL * N_UL);
    Psi_UL{mMTC, episode} = Gamma3 / (beta_mMTC_UL * N_UL);

    % [2.2] Calculate the reward according to equation [13][14]
    R_UL{V2X} = 0;
    R_UL{eMBB} = 0;
    R_UL{mMTC} = 0;
    if Psi_UL{V2X, episode} <= 1
        R_UL{V2X} = exp(Psi_UL{V2X, episode});
    else
        R_UL{V2X} = 1/exp(Psi_UL{V2X, episode}); % TODO: 1/Psi_UL{V2X} or 1/exp(Psi_UL{V2X})?
    end
    if Psi_UL{eMBB, episode} <= 1
        R_UL{eMBB} = exp(Psi_UL{eMBB, episode});
    else
        R_UL{eMBB} = 1/exp(Psi_UL{eMBB, episode}); % TODO: 1/Psi_UL{eMBB} or 1/exp(Psi_UL{eMBB})?
    end
    if Psi_UL{mMTC, episode} <= 1
        R_UL{mMTC} = exp(Psi_UL{mMTC, episode});
    else
        R_UL{mMTC} = 1/exp(Psi_UL{mMTC, episode}); % TODO: 1/Psi_UL{mMTC} or 1/exp(Psi_UL{mMTC})?
    end      
    R_TOT_UL = (R_UL{V2X} * R_UL{eMBB} * R_UL{mMTC}) ^ (1/3);    

    % [3] Q-learning table update according to equation [15]
    if episode > 1
        for ri = 1 : Ar
            for xi = 1 : Ax
                if (ri == ar && xi == ax)
                    Q_UL(ar, ax, episode) = (1 - alpha) * Q_UL(ar, ax, episode-1) + alpha * R_TOT_UL;
                else
                    Q_UL(ri, xi, episode) = Q_UL(ri, xi, episode-1);
                end
            end
        end
    else
        Q_UL(ar, ax, episode) = alpha * R_TOT_UL;
    end
    fprintf(1, 'episode=%04d: ar=%02d, ax=%02d, Gamma1{UL}=%3.2f, Gamma1{SL}=%3.2f, Gamma2=%3.2f, Gamma3=%3.2f, R_UL(V2X)=%3.2f, R_UL(eMBB)=%3.2f, R_UL(mMTC)=%3.2f, R_TOT_UL=%3.2f, beta_V2X=%3.2f, beta_eMBB=%3.2f, beta_mMTC=%3.2f\n', ...
        episode, ar, ax, Gamma1{UL}, Gamma1{SL}, Gamma2, Gamma3, R_UL{V2X}, R_UL{eMBB}, R_UL{mMTC}, R_TOT_UL, beta_V2X_UL, beta_eMBB_UL, beta_mMTC_UL);
end

% To plot the changes of the Q-table during the RL algorithm. 
% For simplicity, only 10 instance of the Q-table during the RL are plotted. 
nr_line_to_plot = 1;
episode_to_check = nr_episode : -nr_episode / nr_line_to_plot : 1;
figure(1);
hold on;
for episode = (episode_to_check)
    mesh(Q_UL(:, :, episode), 'DisplayName', sprintf('episode=%d', episode));     
end
xlabel('actions ax(eMBB, V2X)')
ylabel('actions ar(mMTC)');
zlabel('Q-values');
legend;
grid on;
hold off;

Q_final = Q_UL(:, :, nr_episode);
[~, I] = max(Q_final, [], 'all');
ar_sel = mod((I - 1), Ar) + 1;
ax_sel = ceil(I / 4);
alpha_sel = heuristic_tuning_f(ar_sel, ax_sel, actions_r, actions_x, Psi_UL, beta_UL(:, ar, ax));
alpha_sel_o = alpha_sel;
end




