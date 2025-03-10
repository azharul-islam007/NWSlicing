function [Thr_o] = mapping_SINR2Speff_f(SINR_dB_i, DL_UL_i)
%%
% Mapping from SINR (dB) into the throughput (spectrum efficiency) according
% to the following reference
% Reference: https://www.etsi.org/deliver/etsi_TR/136900_136999/136942/08.02.00_60/tr_136942v080200p.pdf
% Annex A: Link level performance model

SINR_dB = SINR_dB_i;
DL_UL = DL_UL_i;

SINR_dB_min = -10;  % Refer to Table A.1 in the reference
Thr_max = 8.8; % bps/Hz. The maximum spectral efficiency is 8.8bps/Hz as required in the Table in the IEEE Access paper
if strcmpi(DL_UL, 'dl')
    alpha = 0.6; % Refer to Table A.1 in the reference, attenuation factor, representing implementation losses
elseif strcmpi(DL_UL, 'ul')
    alpha = 0.4;
end

SINR_max = 2 ^ (Thr_max) - 1;   % The SNIR at which max throughput is reached according to Shannon bound S(SINR)= log2(1+SNIR) bps/Hz
SINR_dB_max = 10 * log10(SINR_max); % transformed from dB value to the linear value

if SINR_dB < SINR_dB_min        % Refer to Annex A on how to approximate the throughput over a channel with a given SNR when using link adaptation.
    Thr = 0;
elseif SINR_dB < SINR_dB_max
    SINR = 10 ^ (SINR_dB / 10);
    Thr = alpha * log2(1 + SINR);
else
    Thr = Thr_max;
end

Thr_o = Thr;

end