function [alpha_sel_o] = slicing_RL_DL_f(rand_seed_i, nr_episode_i, tao_i, Avg_T_i, V_v2x_i, M_eMBB_i, M_mMTC_i)
%% RAN slicing algorithm, DL

rng(rand_seed_i);   % seed initialization for repeatibility
tao = tao_i;        % The temperature parameter in Table 1 in the paper, refer to section 'SELECTION CRITERION' for its impact on the selection probability
nr_episode = nr_episode_i;  % number of episodes. Every 'action' is evaluated in one episode and a new 'action' will be taken after the completition of an episode.
B = 30e3 * 12;              % bandwidth of a PRB, equal to 30kHz * 12 subcarriers
alpha = 0.1;                % Table-1, learning rate
Sm = 300;                   % Message size, byptes
N_DL = 200;                 % Number of PRBs in DL
Nr_S = 3;                   % 2 service types, {'V2X', 'eMBB'}
Ar = 4;                     % Number of action for slicing of mMTC
Ax = 20;                    % Number of actions for slicing allocation between eMBB and V2X
Fd = 0.5e-3;                % TTI duration, should be 0.5ms because \Delta_f (subcarrier interval) is 30kHz
beta_DL = zeros(Nr_S, Ar, Ax);  % Actions of the RL algorithm

V2X = 1; 
eMBB = 2;
mMTC = 3;
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

Avg_T = Avg_T_i;            % Averaging window before taking the next action
R_DL = cell(Nr_S, 1);       % reward after taking an action
C = 4;                      % number of cluster

% The off-line network model for the V2X serive, which should simulate the
% characteristics of the real network model
% TODO: to integrate the off-line network model with the (simulated)
% on-line network model
V_v2x = V_v2x_i;
pd_v2x_SINR_dB = cell(C, 1);
pd_v2x_Nr_packet = cell(C, 1);
for j = 1 : C
    Vj = V_v2x(j);
    pd_v2x_SINR_dB{j} = cell(Vj, 1);
    pd_v2x_Nr_packet{j} = cell(Vj, 1);
    for i = 1 : Vj
        % off-line network modeling of the V2X SINR. Here is a uniform
        % distribution with a lower bound -10dB and an uppoer bound = 40dB.
        pd_v2x_SINR_dB{j}{i} = makedist('Uniform', 'Lower', -10, 'Upper', 40); 
        % off-line network modeling of the V2X packet arrival rate. Here is
        % a binomial distribution.
        pd_v2x_Nr_packet{j}{i} = makedist('Binomial', 'N', 2, 'p', 0.5);
    end
end

% The off-line network model for the eMBB serive, which should simulate the
% characteristics of the real network model
% TODO: to integrate the off-line network model with the (simulated)
% on-line network model
M_eMBB = M_eMBB_i;
pd_nr_eMBB_session = cell(M_eMBB, 1);
pd_eMBB_SINR_dB = cell(M_eMBB, 1);
for m = 1 : M_eMBB
    % off-line network modeling of the eMBB SINR. Here is a uniform
    % distribution with a lower bound -10dB and an uppoer bound = 40dB.
    pd_eMBB_SINR_dB{m} = makedist('Uniform', 'Lower', -10, 'Upper', 40); 
    % off-line network modeling of the eMBB session arrival rate. Here is a
    % binomial distribution.
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
        pd_nr_eMBB_session, ...
        pd_eMBB_SINR_dB);

    % [1.3] For mMTC serice, to get the number of sporadic transmissions,
    % the packet size and the SINR for the mMTC UEs during the
    % averaging window period from the off-line network model
    [tx_mMTC_DL, ps_mMTC_DL, Speff_mMTC_DL] = offline_nw_mMTC_DL_f(...
        Avg_T, ...
        M_mMTC, ...
        pd_mMTC_SINR_dB, ...
        pd_packet_size_mMTC, ...
        pd_tx_mMTC);

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
            nominator = nominator + tx{t}{m} * ps{t}{m} / Speff_mMTC{t}{m};
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
    fprintf(1, 'episode=%04d: ar=%02d, ax=%02d, Gamma1=%3.2f, Gamma2=%3.2f, Gamma3=%3.2f, R_DL(V2X)=%3.2f, R_DL(eMBB)=%3.2f, R_DL(mMTC)=%3.2f, R_TOT_DL=%3.2f, beta_V2X=%3.2f, beta_eMBB=%3.2f, beta_mMTC=%3.2f\n', ...
        episode, ar, ax, Gamma1, Gamma2, Gamma3, R_DL{V2X}, R_DL{eMBB}, R_DL{mMTC}, R_TOT_DL, beta_V2X_DL, beta_eMBB_DL, beta_mMTC_DL);
end

% To plot the changes of the Q-table during the RL algorithm. 
% For simplicity, only 10 instance of the Q-table during the RL are plotted. 
nr_line_to_plot = 1;
episode_to_check = nr_episode : -nr_episode / nr_line_to_plot : 1;
figure(1);
hold on;
for episode = (episode_to_check)
    mesh(Q_DL(:, :, episode), 'DisplayName', sprintf('episode=%d', episode));   
end
xlabel('actions ax(eMBB, V2X')
ylabel('actions ar(mMTC)');
zlabel('Q-values');
legend;
grid on;
hold off;

% Select the final action (r, k) by selecting the one with the maximum Q
% value
Q_final = Q_DL(:, :, nr_episode);
[~, I] = max(Q_final, [], 'all');
ar_sel = mod((I - 1), Ar) + 1;
ax_sel = ceil(I / 4);
alpha_sel = heuristic_tuning_f(ar_sel, ax_sel, actions_r, actions_x, Psi_DL, beta_DL(:, ar, ax));
alpha_sel_o = alpha_sel;
end




