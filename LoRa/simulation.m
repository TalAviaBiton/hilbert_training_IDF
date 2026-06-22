close all; clear; clc;

SF = [ 7 8 9 10 11 12]; 
BW = [0.5 0.75 1 1.25 1.5] ; % MHz
Pt = 1; % dBM
Fs = BW .* 2; % MHz
SNR = 0; % dB

% Coherece = 1;
% df = 0;
% code_rate = 1;
% symbols_in_preamble = 8;
% sync_key = 5;
% % varargin_Tx = [code_rate, symbols_in_preamble, sync_key];
% Preamble_Symbol_number = symbols_in_preamble;
% % varargin_Rx = [SNR, Preamble_Symbol_number];
% message = [0 0 0 0 1 1 1 1 0 0 0 1 1 1 0 0 1 1 0 1];

t_RS = zeros(length(BW), length(Fs), length(SF));
t_TOA = zeros(length(BW), length(Fs), length(SF));
t_range = zeros(length(BW), length(Fs), length(SF));

for i = 1:length(BW)
    for j = 1:length(SF)
% למה זה דרוש?
        % [signal_mod] = LoRa_Tx(message,BW(i) ,SF(j),Pt,Fs(i),df,code_rate, symbols_in_preamble, sync_key);
        % recived_message = LoRa_Rx(signal_mod,BW(i),SF(j),Coherece,Fs(i),df,SNR, Preamble_Symbol_number);

        C = 3 * 10 ^ 8; % m/s
        lambda = C ./ Fs(i); % m
        N = -174 + BW(i); % dBm
        Pr = N + SNR; % dBm
        Gr = 0; % dBi
        Gt = 0; % dBi
        M = 2 .^ SF(j);

        range = (10 .^ ((Pr - Pt - Gr - Gt) / 20)) * 4 * pi / lambda; % friis, m
        TOA = M.' ./ BW(i); % sec
        Rs = ((BW(i)).' ./ M).'; % per sec

        t_RS(i,i,j) = Rs;
        t_TOA(i,i,j) = TOA;
        t_range(i,i,j) = range;

    end
end


