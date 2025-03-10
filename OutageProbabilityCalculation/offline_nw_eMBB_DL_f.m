function [Rb_o, Speff_o] = offline_nw_eMBB_DL_f(Avg_T_i, ...
    M_eMBB_i, ...
    pd_nr_session_i, ...
    pd_SINR_dB_i ...
    )
% To simulate the eMBB session generation, SINR according to the offline
% network model
% Input: 
%  M_eMBB_i: the number of eMBB UEs in every slot during the averaging
%  window
%  pd_nr_session_i: the random generator for the eMBB session in the
%  offline network model
%  pd_SINR_dB_i: the random generator for the SINR in the offline network
%  model
% Output:
%  Rb_o: the bit rates of every eMBB UE in every slot in the averaging
%  window
%  Speff_o: the spectrum efficiency (i.e. throughput) of every eMBB UE in
%  every slot in the averaging window

Avg_T = Avg_T_i;
M_eMBB = M_eMBB_i;
pd_SINR_dB = pd_SINR_dB_i;
pd_nr_session = pd_nr_session_i;
Rb = cell(Avg_T, 1);
Speff = cell(Avg_T, 1);
for t = 1 : Avg_T
    Rb{t} = cell(M_eMBB, 1);
    Speff{t} = cell(M_eMBB, 1);
    for m = 1 : M_eMBB
        Rb{t}{m} = 1e6 * pd_nr_session{m}.random(); % one eMBB session requires 1Mbps throughput
        SINR_dB = pd_SINR_dB{m}.random();
        Speff{t}{m} = mapping_SINR2Speff_f(SINR_dB, 'dl');
    end
end

Rb_o = Rb;
Speff_o = Speff;
end