function [m_DL_o, ...
    Speff_DL_o] = offline_nw_v2x_DL_f(AvgT_i, ...
    C_i, ...
    V_i, ...
    pd_SINR_dB_i, ...
    pd_nr_packet_i)

% To simulate the V2X packets generation, SINR according to the offline
% network model
% Input: 
%  C_i: the number of clusters 
%  V_i: the number of vehicular UE in every cluster
%  pd_SINR_dB_i: the random generator for the SINR in the offline network
%  model
%  pd_nr_packet_i: the random generator for the V2X packets in the
%  offline network model
% Output:
%  m_DL_o: the number of transmitted packets by the vehicular UE in every
%  slot in the averaging window 
%  Speff_DL_o: the spectrum efficiency (i.e.
%  throughput) of every vehicular UE in every slot in the averaging window

AvgT = AvgT_i;
C = C_i;
V = V_i;
pd_SINR_dB = pd_SINR_dB_i;
pd_nr_packet = pd_nr_packet_i;
m_DL = cell(AvgT_i, 1);
Speff_DL = cell(AvgT_i, 1);
 
for t = 1 : AvgT
    m_DL{t} = cell(C, 1);
    Speff_DL{t} = cell(C, 1);
    for j = 1 : C
        Vj = V(j);
        m_DL{t}{j} = cell(Vj, 1);
        Speff_DL{t}{j} = cell(Vj, 1);
        for i = 1 : V(j)
            m_DL{t}{j}{i} = pd_nr_packet{j}{i}.random();
            SINR_dB = pd_SINR_dB{j}{i}.random();
            Speff_DL{t}{j}{i} = mapping_SINR2Speff_f(SINR_dB, 'dl');
        end
    end
end

m_DL_o = m_DL;
Speff_DL_o = Speff_DL;

end