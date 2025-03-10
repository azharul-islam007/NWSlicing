function [alpha_o] = heuristic_tuning_f(...
    ar_sel_i, ax_sel_i, ...
    actions_r_i, actions_x_i, ...
    Psi_i, beta_i)

%% LOW-COMPLEXITY HEURISTICS SOLUTION
% Input:
% ar_sel_i: the selection of mMTC
% ax_sel_i: the selection of eMBB and V2X;
% actions_r_i: the list of all the actions selected for mMTC during the RL training
% actions_x_i: the list of all the actions selected for eMBB and V2X during the RL
% training
% Psi_i: the list of all the normalized resource utilization during the RL training
% beta_i: the action sets (refer to the report for its definition)

ar_sel = ar_sel_i;
ax_sel = ax_sel_i;
actions_r = actions_r_i;
actions_x = actions_x_i;
Psi = Psi_i;
Nr_s = size(Psi, 1); 
V2X = 1;
eMBB = 2;
mMTC = 3;

nr_episode = size(actions_r, 1);    
% To estimate the average of the Psi (PRB utilization for each service) for
% the selected action (ar, ax) using the records acquired from the RL
% training
Psi_sel_sum = zeros(Nr_s, 1); 
count_sel = 0;
for episode = 1 : nr_episode
    if actions_r(episode) == ar_sel && actions_x(episode) == ax_sel
        Psi_sel_sum(V2X) = Psi_sel_sum(V2X) + Psi{V2X, episode};
        Psi_sel_sum(eMBB) = Psi_sel_sum(eMBB) + Psi{eMBB, episode};
        Psi_sel_sum(mMTC) = Psi_sel_sum(mMTC) + Psi{mMTC, episode};
        count_sel = count_sel + 1;
    end
end
Psi_sel_avg = Psi_sel_sum / count_sel;

% fine tuning the \beta(s,x) into \alpha(s,x) according the algorithm in the report
Delta_C = zeros(Nr_s, 1);
beta = beta_i;
omega = 0.85; % fixed according to the paper
W = 0;
for s = [V2X eMBB mMTC]
    if Psi_sel_avg(s) < 1
        Delta_C(s) = (1 - Psi_sel_avg(s)) * omega;
    else
        W = W + Psi_sel_avg(s);
    end
end
Delta_C_sum = sum(Delta_C);

if (W > 0  && Delta_C_sum > 0)
    alpha = zeros(Nr_s, 1);
    for s = [V2X eMBB mMTC]
        if Psi_sel_avg(s) < 1
            alpha(s) = beta(s) - Delta_C(s); % reduce a fraction  of extra capacity
        else
            alpha(s) = beta(s) + Delta_C_sum * Psi_sel_avg(s) / W; % add a fraction of extra capacity
        end
    end
else
    alpha = beta;
end

alpha_o = alpha;
end