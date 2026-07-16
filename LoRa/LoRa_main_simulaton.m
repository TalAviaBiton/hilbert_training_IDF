%% define SF & BW and get a BET curve and a spectrogram of the transmission as well as the range TOA and Rs 
% close all;
clear; clc;

%% params definitions
message = 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII';
EbNo_vec = -10:2:10; % dB
SF =  7 ; 
BW = 1.5 * 1e6; % Hz
Pt = -40; % dBm 
Fs = BW .* 2; % Hz
C = 3 * 10 ^ 8; % m/s
Gr = 0; % dBi
Gt = 0; % dBi
Coherence = 1; 
df = 0;
code_rate = 1;
symbols_in_preamble = 8;
sync_key = 5;
Preamble_Symbol_number = symbols_in_preamble;

%% sending the signal
lambda = C / Fs; % m
N_noise = -174 + 10*log10(BW); % dBm/Hz noise figure
[signal_mod] = LoRa_Tx(message, BW, SF, Pt, Fs, df, code_rate, symbols_in_preamble, sync_key);

%% displaying in spectrogram
figure('Name', 'Spectrogram of Received LoRa Signal');
spectrogram(signal_mod, 256, 250, 256, Fs, 'yaxis');
title('Spectrogram of Modulated Signal');

%% Pre-calculating to save time
msg_tx_char = convertStringsToChars(message);
msg_tx_bits = reshape(dec2bin(msg_tx_char, 8).', 1, []).'- '0';
ber_results = zeros(size(EbNo_vec));

snr_db_vec = EbNo_vec + 10*log10(code_rate*SF/(2^SF));

%% displaying BER curve
for k = 1:length(EbNo_vec)
    % receive signal using pre-calculated values
    recived_message = LoRa_Rx(signal_mod, BW, SF, Coherence, Fs, df, snr_db_vec(k), Preamble_Symbol_number);
    msg_rx_char = convertStringsToChars(recived_message);
    
    % calculate BER
    if length(msg_rx_char) == length(msg_tx_char)
        msg_rx_bits = reshape(dec2bin(msg_rx_char, 8).', 1, []).'- '0';
        [~, ber_results(k)] = biterr(msg_tx_bits, msg_rx_bits);
    else
        ber_results(k) = 0.5; 
    end
end

% display
figure('Name', ['BER Curve - SF' num2str(SF) ' BW' num2str(BW/1e6) 'MHz']);
semilogy(EbNo_vec, ber_results, '-o', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title(['BER Performance - SF = ' num2str(SF) ', BW = ' num2str(BW/1e6) ' MHz']);

%% calculate results
M = 2 ^ SF;
Pr_vec = N_noise + snr_db_vec; 

TOA = M / BW; % sec
RS = BW / M;  % symbols per sec
range = ((10 .^ ((Pr_vec - Pt - Gr - Gt) / 20)) * 4 * pi / lambda).'; % meters
