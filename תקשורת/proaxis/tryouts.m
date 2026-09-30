close all; clc; clear;

%% Set Parameters
M = 2;                  % Modulation order (e.g., 2 for BPSK, 4 for QPSK)
pulse_name = 'rect';    % 'rc', 'rrc', or 'rect'

Rs = 1e5;               % Symbol rate (100 kBaud)
fc = 4e5;               % Carrier frequency (400 kHz)
fs = 4 * fc;            % Sampling frequency (must be > 2*fc for Passband)
sps = fs / Rs;          % Samples per symbol

num_symbols = 100;
num_bits = num_symbols * log2(M); 
data_bits = randi([0 1], 1, num_bits); % Generate random bits

beta = 0.5;             % Roll-off factor for RC/RRC
span = 6;               % Filter span in symbols

%% Run MPSK Transmitter
[MPSK_passband, MPSK_baseband, t] = MPSK(data_bits, M, pulse_name, beta, span, sps, fc, fs);

%% Plot Constellation
% Scatterplot of the baseband signal (downsampled to symbol rate)
scatterplot(MPSK_baseband(1:sps:end));
title(sprintf('%d-PSK Baseband Constellation', M));


%% Main Functions
function [BP_signal, BB_signal, t] = MPSK(data_bits, M, pulse_name, beta, span, sps, fc, fs)
    % 1. Bits to Symbols
    % Make sure bits are shaped correctly for bi2de
    data_bits_matrix = reshape(data_bits.', log2(M), []).';
    data_symbols = bi2de(data_bits_matrix, 'left-msb');
    
    % Get complex constellation points
    constellations = get_constellations(M);
    modulated_data = constellations(data_symbols + 1);
    
    % 2. Upsample the data (insert zeros between symbols)
    upsampled_data = upsample(modulated_data, sps);
    
    % 3. Pulse Shaping
    pulse = get_pulse(pulse_name, beta, span, sps);
    % Use 'same' to keep signal length aligned
    BB_signal = conv(upsampled_data, pulse, 'same'); 
    
    % 4. Baseband to Passband
    [BP_signal, t] = BB2BP(BB_signal, fc, fs);
end

function constellations = get_constellations(M)
    % Correct check for power of 2
    if mod(log2(M), 1) ~= 0
        error("M must be a power of 2");
    end
    
    shift = 2 * pi / M;
    phases = zeros(M, 1);
    
    % Optional: Rotate QPSK so it sits on diagonals like standard LTE/Wi-Fi
    if M == 4
        phases(1) = pi / 4; 
    end
    
    for i = 2:M
        phases(i) = phases(i - 1) + shift;
    end
    
    % Convert angles to complex exponential (I + jQ)
    constellations = exp(1j * phases);
end

function pulse = get_pulse(pulse_name, beta, span, sps)
    if strcmp(pulse_name, 'rc')
        pulse = rcosdesign(beta, span, sps, 'normal');
    elseif strcmp(pulse_name, 'rrc')
        pulse = rcosdesign(beta, span, sps, 'sqrt');
    elseif strcmp(pulse_name, 'rect')
        % A simple rectangular pulse normalized by energy
        pulse = ones(sps, 1) / sqrt(sps);
    else
        error("pulse name must be rc / rrc / rect");
    end
end

function [BP_signal, t] = BB2BP(BB_signal, fc, fs)
    % Generate discrete time vector perfectly matched to the signal length
    t = (0:length(BB_signal)-1).' / fs;
    
    % CRITICAL FIX: Passband is the Real part of the mixed signal
    % Re{ Baseband * exp(j * 2*pi * fc * t) }
    BP_signal = real(BB_signal .* exp(1j * 2 * pi * fc * t));
    
    % Plot the first 200 samples to zoom in on the waveform
    num_samples_to_plot = min(200, length(t));
    t_plot = t(1:num_samples_to_plot);
    
    figure('Name', 'Waveforms (Zoomed In)');
    subplot(2,1,1);
    plot(t_plot, real(BB_signal(1:num_samples_to_plot)), 'LineWidth', 1.5);
    title('Baseband Signal (Real Part / In-Phase)');
    xlabel('Time (s)'); grid on;
    
    subplot(2,1,2);
    plot(t_plot, BP_signal(1:num_samples_to_plot), 'r', 'LineWidth', 1.5);
    title('Passband Signal (Upconverted)');
    xlabel('Time (s)'); grid on;
end