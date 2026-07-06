%% define SF & BW and get a BET curve and a spectrogram of the transmission as well as the range TOA and Rs 
% close all;
clear; clc;

%% params definitions
%LORAPHY LoRa physical layer implementation
%%% Example %%%
rf_freq = 470e6;
sf = 12;
bw = 1.25e6;
fs = 1e6;
phy = LoRaPHY(rf_freq, sf, bw, fs);
phy.has_header = 1; % explicit header mode

message = [0 0 0 0 1 1 1 1 0 0 0 1 1 1 0 0 1 1 0 1];
symbols = phy.encode(message);

%% displaying in spectrogram
    phy.spec(symbols, fs, bw, sf);

%% Pre-calculating to save time
EbNo_vec = -10:2:10; % dB
ber_results = zeros(size(EbNo_vec));
N_noise = -174 + 10*log10(bw); % dBm/Hz noise figure
snr_db_vec = EbNo_vec + 10*log10(sf/(2^sf));

%% displaying BER curve
for k = 1:length(EbNo_vec)

    noise = randi(0:1, 1,length(data));
    [data, ~] = phy.decode(symbols + noise);

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
