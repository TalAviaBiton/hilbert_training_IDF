%% define SF & BW and get a BET curve and a spectrogram of the transmission as well as the range TOA and Rs 
close all;
clear; clc;

%% params definitions
%LORAPHY LoRa physical layer implementation
%%% Example %%%
rf_freq = 470e6;
sf = 7;
bw = 125e3;
fs = 1e6;
phy = LoRaPHY(rf_freq, sf, bw, fs);
phy.has_header = 1; % explicit header mode

fs = 10000;             % תדר דגימה (Hz)
duration_ms = 50;       % משך האות במילי-שניות
t = 0:1/fs:(duration_ms/1000 - 1/fs); 

f0 = 100;               % תדר התחלתי
f1 = 2000;              % תדר סופי

y = chirp(t, f0, duration_ms/1000, f1, 'linear');

bits = y > 0;

% t = 0:1/20e3:2;
% chirpy = chirp(t,20,100,250);
% pspectrum(chirpy,1e3,"spectrogram",TimeResolution=0.1, ...
%     OverlapPercent=99,Leakage=0.85)

% inst_freq = instfreq(chirpy, fs);
% mid_freq = (max(inst_freq) + min(inst_freq)) / 2;
% message = inst_freq > mid_freq; 
symbols = phy.encode(bits);

%% displaying in spectrogram
    phy.spec(symbols, fs, bw, sf);

%% Pre-calculating to save time
EbNo_vec = -10:2:10; % dB
ber_results = zeros(size(EbNo_vec));
N_noise = -174 + 10*log10(bw); % dBm/Hz noise figure
snr_db_vec = EbNo_vec + 10*log10(sf/(2^sf));

%% displaying BER curve
for k = 1:length(EbNo_vec)

    % add noise
    SNR_lin = 10.^(snr_db_vec/10);
    P_sig = mean(abs(symbols).^2);
    P_noise = P_sig ./ SNR_lin;
    noisy_symbols = symbols + sqrt(P_noise/2) .* randn(size(symbols));
    
    [data, ~] = phy.decode(noisy_symbols);

    % calculate BER
    if length(data) == length(message)
        [~, ber_results(k)] = biterr(data, message);
    else
        ber_results(k) = 0.5; 
    end
end

% display
figure('Name', ['BER Curve - SF' num2str(sf) ' BW' num2str(bw/1e6) 'MHz']);
semilogy(EbNo_vec, ber_results, '-o', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title(['BER Performance - SF = ' num2str(sf) ', BW = ' num2str(bw/1e6) ' MHz']);

%% calculate results
M = 2 ^ sf;
Pr_vec = N_noise + snr_db_vec; 

TOA = M / bw; % sec
RS = bw / M;  % symbols per sec

C = 3 * 10 ^ 8; % m/s
lambda = C / fs; % m

range = ((10 .^ (Pr_vec  / 20)) * 4 * pi / lambda).'; % meters
