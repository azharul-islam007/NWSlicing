function [m_UL_o, ...
    m_SL_o, ...
    Speff_UL_o, ...
    Speff_SL_o] = offline_nw_v2x_UL_f(AvgT_i, ...
    C_i, ...
    Vj_UL_i, ...
    pd_UL_SINR_dB_i, ...
    pd_UL_Nr_packet_i, ...
    Vj_SL_i, ...
    pd_SL_SINR_i, ...
    pd_SL_Nr_packet_i)

% To simulate the V2X packets generation, SINR according to the offline
% network model
% Input: 
%  C_i: the number of clusters 
%  Vj_UL_i: the number of vehicular UE in every cluster under cellular mode
%  pd_UL_SINR_dB_i: the random generator for the SINR in the offline network
%  model for vehicular UE in cellular mode
%  pd_UL_SINR_dB_i: the random generator for the V2X packets in the
%  offline network model for vehicular UE in cellular mode
%  Vj_SL_i: the number of vehicular UE in every cluster under sidelink mode
%  pd_SL_SINR_dB_i: the random generator for the SINR in the offline network
%  model for vehicular UE in sidelink mode
%  pd_SL_SINR_dB_i: the random generator for the V2X packets in the
%  offline network model for vehicular UE in sidelink mode
% Output:
%  m_UL_o: the number of transmitted packets by the vehicular UE in
%  cellular mode in every slot in the averaging window 
%  m_SL_o: the number of transmitted packets by the vehicular UE in
%  sidelink mode in every slot in the averaging window 
% Speff_UL_o: the spectrum efficiency (i.e. throughput) of every vehicular
%  UE in cellular mode in every slot in the averaging window
% Speff_SL_o: the spectrum efficiency (i.e. throughput) of every vehicular
%  UE in sidelink mode in every slot in the averaging window

AvgT = AvgT_i;
C = C_i;
Vj_UL = Vj_UL_i;
pd_UL_SINR_dB = pd_UL_SINR_dB_i;
pd_UL_Nr_packet = pd_UL_Nr_packet_i;
m_UL = cell(AvgT_i, 1);
Speff_v2x_UL = cell(AvgT_i, 1);

% cellular mode
for t = 1 : AvgT
    m_UL{t} = cell(C, 1);
    Speff_v2x_UL{t} = cell(C, 1);
    for j = 1 : C
        Vj = Vj_UL(j);
        m_UL{t}{j} = cell(Vj, 1);
        Speff_v2x_UL{t}{j} = cell(Vj, 1);
        for i = 1 : Vj_UL(j)
            m_UL{t}{j}{i} = pd_UL_Nr_packet{j}{i}.random();
            SINR_dB = pd_UL_SINR_dB{j}{i}.random();
            Speff_v2x_UL{t}{j}{i} = mapping_SINR2Speff_f(SINR_dB, 'ul');
        end
    end
end

% sidelink mode
Vj_SL = Vj_SL_i;
pd_SL_SINR = pd_SL_SINR_i;
pd_SL_Nr_packet = pd_SL_Nr_packet_i;
m_SL = cell(AvgT_i, 1);
Speff_v2x_SL = cell(AvgT_i, 1);

for t = 1 : AvgT
    m_SL{t} = cell(C, 1);
    Speff_v2x_SL{t} = cell(C, 1);
    for j = 1 : C
        Vj = Vj_SL(j);
        m_SL{t}{j} = cell(Vj, 1);
        Speff_v2x_SL{t}{j} = cell(Vj, 1);
        for i = 1 : Vj_SL(j)
            m_SL{t}{j}{i} = pd_SL_Nr_packet{j}{i}.random();
            SINR_dB = pd_SL_SINR{j}{i}.random();
            Speff_v2x_SL{t}{j}{i} = mapping_SINR2Speff_f(SINR_dB, 'ul');
        end
    end
end

m_UL_o = m_UL;
m_SL_o = m_SL;
Speff_UL_o = Speff_v2x_UL;
Speff_SL_o = Speff_v2x_SL;


end