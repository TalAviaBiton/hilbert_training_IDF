% close all; clear; clc;

SF = [7 8 9 10 11 12]; 
BW = [0.5 0.75 1 1.25 1.5] *1e6; % Hz
Pt = 1; % dBM
Fs = BW .* 2; % MHz
SNR = 100; % dB
C = 3 * 10 ^ 8; % m/s
Gr = 0; % dBi
Gt = 0; % dBi
Coherece = 1;
df = 0;
code_rate = 1;
symbols_in_preamble = 8;
sync_key = 5;
% varargin_Tx = [code_rate, symbols_in_preamble, sync_key];
Preamble_Symbol_number = symbols_in_preamble;
% varargin_Rx = [SNR, Preamble_Symbol_number];
message = [0 0 0 0 1 1 1 1 0 0 0 1 1 1 0 0 1 1 0 1];

l_BW = length(BW);
l_SF = length(SF);

t_RS = zeros(l_BW, l_SF);
t_TOA = zeros(l_BW, l_SF);
t_range = zeros(1, l_BW);

for i = 1:l_BW

    lambda = C ./ Fs(i); % m
    N = -174 + BW(i); % dBm
    Pr = N + SNR; % dBm
    t_range(i) = (10 .^ ((Pr - Pt - Gr - Gt) / 20)) * 4 * pi / lambda; % friis, m

    for j = 1:l_SF

        [signal_mod] = LoRa_Tx(message,BW(i) ,SF(j),Pt,Fs(i),df,code_rate, symbols_in_preamble, sync_key);
        recived_message = LoRa_Rx(signal_mod,BW(i),SF(j),Coherece,Fs(i),df,SNR, Preamble_Symbol_number);

        M = 2 .^ SF(j);

        t_TOA(i,j) = M.' ./ BW(i); % sec
        t_RS(i,j) = ((BW(i)).' ./ M).'; % per sec

    end
end


