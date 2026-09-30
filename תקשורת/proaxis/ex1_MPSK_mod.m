close all; clc; clear;

% set params:
M = 8;
pulse_name = 'rect';

Rs = 1e5;               % Symbol rate (100 kBaud)
fc = 1e6;
fs = 4 * fc;
sps = fs / Rs;          % Samples per symbol

num_symbols = 100;
num_bits = num_symbols * log2(M); 
data_bits = randi([0 1], 1, num_bits); % Generate random bits

beta = 0.5;             % Roll-off factor for RC/RRC
span = 6;               % Filter span in symbols
bw = beta;

[MPSK_passband, MPSK_baseband, t] = MPSK(data_bits, M, pulse_name, beta, span, sps, fc, fs, bw);

scatterplot(MPSK_baseband(1:sps:end));
title(sprintf('%d-PSK Baseband Constellation', M));

function [BP_signal, BB_signal, t] = MPSK(data_bits, M, pulse_name, beta, span, sps, fc, fs, bw)
data_bits_matrix = reshape(data_bits.', log2(M), []).';
data_symbols = bi2de(data_bits_matrix, 'left-msb');
constellations = get_constellations(M);
modulated_data = constellations(data_symbols + 1);
upsampled_data = upsample(modulated_data, sps);
t_pulse = (-span*sps/2 : span*sps/2).' / fs;
pulse = get_pulse(pulse_name, beta, span, sps, t_pulse, fc, bw);
BB_signal = conv(upsampled_data, pulse, 'same');
[BP_signal, t] = BB2BP(BB_signal, fc, fs);
end

function constellations = get_constellations(M)

 if mod(log2(M), 1) ~= 0
        error("M must be a power of 2");
 end
    
shift = 2 * pi / M;
phases = zeros(M, 1);

if M == 4
        phases(1) = pi / 4; 
end

for i = 2:M
    phases(i, 1) = phases(i - 1) + shift;
end

constellations = exp(1j * phases);

end

function pulse = get_pulse(pulse_name, beta, span, sps, t, fc, bw)

if strcmp(pulse_name, 'rc')
    pulse = rcosdesign(beta, span, sps);
elseif strcmp(pulse_name, 'rrc')
    pulse = rcosdesign(beta, span, sps, "sqrt");
elseif strcmp(pulse_name, 'gaussian')
    [~,~,pulse] = gauspuls(t, fc, bw);
elseif strcmp(pulse_name, 'rect')
    pulse = ones(sps, 1) / sqrt(sps);
else
    error("pulse name must be rc/rrc/rect/gaussian");
end

end

function [BP_signal, t] = BB2BP(BB_signal, fc, fs)

    t = (0:length(BB_signal)-1).' / fs;
    BP_signal = real(BB_signal .* exp(2 * 1j * pi * fc * t));

    %plots
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




