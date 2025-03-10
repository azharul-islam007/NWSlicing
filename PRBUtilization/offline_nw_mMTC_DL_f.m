function [tx_DL_o, ...
    ps_DL_o, ...
    Speff_DL_o] = offline_nw_mMTC_DL_f(AvgT_i, ...
    M_mMTC_i, ...
    pd_SINR_dB_i, ...
    pd_log2_packet_size_i, ...
    pd_tx_i)

% To simulate the mMTC session generation, SINR according to the offline
% network model 
% Input:
%  M_mMTC_i: the number of mMTC UEs in every slot during the averaging
%  window
%  pd_SINR_dB_i: the random generator for the SINR in the offline network
%  pd_packet_size_i: the random generator for the mMTC packet size in the
%  offline network model
%  pd_tx_i: the random generator for the transmission probability of the
%  mMTC UEs.
%  model
% Output:
%  tx_DL_o: the transmission indication flag in every slot for every mMTC
%  UE
%  ps_DL_o: the packet size of every mMTC UE in every slot in the averaging
%  window
%  Speff_DL_o: the spectrum efficiency (i.e. throughput) of every mMTC UE in
%  every slot in the averaging window

AvgT = AvgT_i;
M_mMTC = M_mMTC_i;
pd_SINR_dB = pd_SINR_dB_i;
pd_log2_packet_size = pd_log2_packet_size_i;
pd_tx = pd_tx_i;

tx_DL = cell(AvgT_i, 1);
ps_DL = cell(AvgT_i, 1);
Speff_DL = cell(AvgT_i, 1);

for t = 1 : AvgT
    tx_DL{t} = cell(M_mMTC, 1);
    ps_DL{t} = cell(M_mMTC, 1);
    Speff_DL{t} = cell(M_mMTC, 1);
    for m = 1 : M_mMTC
        tx_DL{t}{m} = pd_tx{m}.random();
        ps_DL{t}{m} = 2 ^ (pd_log2_packet_size{m}.random());
        SINR_dB = pd_SINR_dB{m}.random();
        Speff_DL{t}{m} = mapping_SINR2Speff_f(SINR_dB, 'dl');
    end
end

tx_DL_o = tx_DL;
ps_DL_o = ps_DL;
Speff_DL_o = Speff_DL;

end