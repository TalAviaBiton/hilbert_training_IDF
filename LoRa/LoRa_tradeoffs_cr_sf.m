% LoRa CSS PHY Simulation: Configurable SF and CR
clear; clc; close all;

%% 1. System Parameters (Change these to explore tradeoffs!)
SF = 7;                     % Spreading Factor 
CR_den = 6;                 % Code Rate denominator: 5, 6, 7, or 8 (for 4/5, 4/6, 4/7, 4/8)
num_symbols = 500;          % Number of LoRa symbols to simulate
BW = 125e3;                 % Bandwidth (Hz)
SNR_dB = 30;                 % Additive White Gaussian Noise level

% Derived Parameters
M = 2^SF;                   % Chips per symbol (FFT size)
T_sym = M / BW;             % Time duration of one symbol
bit_rate = SF * (4/CR_den) / T_sym; % Effective physical data rate

%% 2. Data Generation & FEC Encoding (Configurable CR)
num_data_blocks = floor((num_symbols * SF) / CR_den);
tx_data = randi([0 1], num_data_blocks, 4); % 4-bit data words

% Generate standard Hamming matrices for LoRa
% Base (7,4) Parity Matrix
P = [1 1 0; 0 1 1; 1 1 1; 1 0 1]; 
G7 = [eye(4) P]; 
% Extend to (8,4) by adding a row-parity bit
G8 = [G7 mod(sum(G7,2), 2)]; 

% Slice the Generator matrix based on selected Code Rate
G = G8(:, 1:CR_den); 

% Encode Data
tx_encoded = mod(tx_data * G, 2);
tx_encoded_stream = tx_encoded.';
tx_encoded_stream = tx_encoded_stream(:);

% Pad with zeros to fill exact number of SF symbols
pad_len = mod(SF - mod(length(tx_encoded_stream), SF), SF);
if pad_len > 0
    tx_encoded_stream = [tx_encoded_stream; zeros(pad_len, 1)];
end

%% 3. Interleaving
num_cols = length(tx_encoded_stream) / SF;
% Interleave: write columns of length SF, read rows
interleaved_bits = reshape(reshape(tx_encoded_stream, SF, num_cols)', [], 1);

%% 4. LoRa CSS Modulation
sym_tx_int = bit2int(reshape(interleaved_bits, SF, []), SF);

n = (0:M-1).';
base_upchirp = exp(1j * pi * ((n.^2) / M - n));
base_downchirp = conj(base_upchirp);

tx_signal = zeros(M * num_symbols, 1);
for k = 1:num_symbols
    shifted_chirp = circshift(base_upchirp, sym_tx_int(k));
    tx_signal((k-1)*M + 1 : k*M) = shifted_chirp;
end

%% 5. Channel: AWGN + Burst/Colored Noise
rx_signal = awgn(tx_signal, SNR_dB, 'measured');

% Inject Burst Noise (corrupts a chunk of adjacent symbols)
burst_start_sym = floor(num_symbols / 4); 
burst_len_sym = 4; % Number of symbols completely wiped out
burst_idx = (burst_start_sym*M) : ((burst_start_sym + burst_len_sym)*M - 1);

% Colored noise (Low-pass filtered complex Gaussian)
[b, a] = butter(4, 0.05); 
colored_noise = filter(b, a, randn(length(burst_idx), 1) + 1j*randn(length(burst_idx), 1));
colored_noise = 10 * colored_noise / std(colored_noise); 
% rx_signal(burst_idx) = rx_signal(burst_idx) + colored_noise;
rx_signal(burst_idx) = colored_noise;

%% 6. Receiver: Demodulation
rx_sym_int = zeros(num_symbols, 1);
for k = 1:num_symbols
    chunk = rx_signal((k-1)*M + 1 : k*M);
    dechirped = chunk .* base_downchirp;
    X = fft(dechirped);
    [~, peak_idx] = max(abs(X));
    rx_sym_int(k) = peak_idx - 1; 
end

%% 7. De-Interleaving & Maximum Likelihood FEC Decoding
rx_interleaved_bits = int2bit(rx_sym_int.', SF);
rx_interleaved_bits = rx_interleaved_bits(:);

% De-interleave
rx_encoded_stream = reshape(reshape(rx_interleaved_bits, num_cols, SF)', [], 1);
burst_errors = rx_interleaved_bits ~= interleaved_bits; % For plotting
scattered_errors = rx_encoded_stream ~= tx_encoded_stream; % For plotting

% Remove padding and reshape back to blocks of size CR_den
rx_encoded_stream = rx_encoded_stream(1:end-pad_len);
rx_encoded = reshape(rx_encoded_stream, CR_den, [])';

% --- Maximum Likelihood Block Decoder ---
% Generate all 16 possible valid codewords for the chosen Code Rate
all_data = dec2bin(0:15) - '0';
valid_codewords = mod(all_data * G, 2);

rx_data = zeros(num_data_blocks, 4);
for i = 1:num_data_blocks
    % Calculate Hamming distance between received block and all 16 valid codewords
    dists = sum(rx_encoded(i,:) ~= valid_codewords, 2);
    [~, min_idx] = min(dists); % Pick the closest match
    rx_data(i,:) = all_data(min_idx, :);
end

% %% 8. Visualization & Tradeoff Summary
% figure('Name', 'LoRa Tradeoffs', 'Position', [100 100 900 600]);
% 
% subplot(3,1,1);
% stem(burst_errors, 'r', 'Marker', '.', 'LineStyle', 'none');
% title('Errors BEFORE De-Interleaving (Burst Noise impact)');
% xlim([burst_start_sym*SF - 50, (burst_start_sym + burst_len_sym)*SF + 50]); grid on;
% 
% subplot(3,1,2);
% stem(scattered_errors, 'b', 'Marker', '.', 'LineStyle', 'none');
% title('Errors AFTER De-Interleaving (Errors are scattered 1 per block)');
% xlim([0, length(scattered_errors)]); grid on;
% 
% subplot(3,1,3);
 tx_bits = tx_data.'; tx_bits = tx_bits(:);
 rx_bits = rx_data.'; rx_bits = rx_bits(:);
 final_errors = tx_bits ~= rx_bits;
% 
% stem(final_errors, 'k', 'Marker', '.', 'LineStyle', 'none');
% title(sprintf('Final Output Errors (CR = 4/%d) | Total: %d', CR_den, sum(final_errors)));
% xlim([0, length(final_errors)]); grid on;

% --- Console Tradeoff Report ---
fprintf('\n=== LORA PHY TRADEOFF REPORT ===\n');
fprintf('Spreading Factor : %d\n', SF);
fprintf('Code Rate        : 4/%d\n', CR_den);
fprintf('Symbol Time      : %.2f milliseconds\n', T_sym * 1000);
fprintf('Data Rate        : %.2f bits/sec\n', bit_rate);
fprintf('--------------------------------\n');
fprintf('Burst Errors     : %d bits\n', sum(burst_errors));
fprintf('Final Errors     : %d bits\n', sum(final_errors));

if sum(final_errors) == 0
    fprintf('FEC Result       : SUCCESS (Redundancy was sufficient to correct errors)\n');
else
    fprintf('FEC Result       : FAILED (Code Rate too weak to correct the scattered errors)\n');
end
fprintf('================================\n');