function [tx_UL_o, ...
    ps_UL_o, ...
    Speff_UL_o] = offline_nw_mMTC_UL_f(AvgT_i, ...
    M_mMTC_i, ...
    pd_SINR_dB_i, ...
    pd_packet_size_i, ...
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
%  tx_UL_o: the transmission indication flag in every slot for every mMTC
%  UE
%  ps_UL_o: the packet size of every mMTC UE in every slot in the averaging
%  window
%  Speff_UL_o: the spectrum efficiency (i.e. throughput) of every mMTC UE in
%  every slot in the averaging window

AvgT = AvgT_i;
M_mMTC = M_mMTC_i;
pd_SINR_dB = pd_SINR_dB_i;
pd_packet_size = pd_packet_size_i;
pd_tx = pd_tx_i;

tx_UL = cell(AvgT_i, 1);
ps_UL = cell(AvgT_i, 1);
Speff_UL = cell(AvgT_i, 1);

for t = 1 : AvgT
    tx_UL{t} = cell(M_mMTC, 1);
    ps_UL{t} = cell(M_mMTC, 1);
    Speff_UL{t} = cell(M_mMTC, 1);
    for m = 1 : M_mMTC
        tx_UL{t}{m} = pd_tx{m}.random();
        ps_UL{t}{m} = pd_packet_size{m}.random();
        SINR_dB = pd_SINR_dB{m}.random();
        Speff_UL{t}{m} = mapping_SINR2Speff_f(SINR_dB, 'ul');
    end
end

tx_UL_o = tx_UL;
ps_UL_o = ps_UL;
Speff_UL_o = Speff_UL;

end