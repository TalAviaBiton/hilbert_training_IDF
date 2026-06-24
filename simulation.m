
close all; clear; clc;

SF = [5 6 7 8 9 10 11 12]; 
BW = [0.5 0.75 1 1.25 1.5] * 1e6; % Hz
Pt = 1; % dBm 
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

message = 'Hello LoRa';
l_BW = length(BW);
l_SF = length(SF);

t_RS = zeros(l_BW, l_SF);
t_TOA = zeros(l_BW, l_SF);
t_range = zeros(l_BW, l_SF);

EbNo_vec = -10:2:10; % dB

for i = 1:l_BW
    lambda = C ./ Fs(i); % m
    N_noise = -174 + 10*log10(BW(i)); % dBm/Hz 

    for j = 1:l_SF
        [signal_mod] = LoRa_Tx(message, BW(i), SF(j), Pt, Fs(i), df, code_rate, symbols_in_preamble, sync_key);

        figure('Name', 'Spectrogram of Received LoRa Signal');
        spectrogram(signal_mod, 256, 250, 256, Fs(i), 'yaxis');
        title('Spectrogram of Last Modulated Signal');

        msg_tx_char = convertStringsToChars(message);
        msg_tx_bits = reshape(dec2bin(msg_tx_char, 8).', 1, []).'- '0';

        ber_results = zeros(size(EbNo_vec));

        for k = 1:length(EbNo_vec)
            snr_db = EbNo_vec(k) + 10*log10(code_rate/4) + 10*log10(SF(j)); 

            Pr = N_noise + snr_db; 
            t_range(i, j) = (10 .^ ((Pr - Pt - Gr - Gt) / 20)) * 4 * pi / lambda;

            recived_message = LoRa_Rx(signal_mod, BW(i), SF(j), Coherence, Fs(i), df, snr_db, Preamble_Symbol_number);

            msg_rx_char = convertStringsToChars(recived_message);

            if length(msg_rx_char) == length(msg_tx_char)
                msg_rx_bits = reshape(dec2bin(msg_rx_char, 8).', 1, []).'- '0';
                [~, ber_results(k)] = biterr(msg_tx_bits, msg_rx_bits);
            else
                ber_results(k) = 0.5; 
            end
        end

        figure('Name', ['BER Curve - SF' num2str(SF(j)) ' BW' num2str(BW(i)/1e6) 'MHz']);
        semilogy(EbNo_vec, ber_results, '-o', 'LineWidth', 2);
        grid on;
        xlabel('E_b/N_0 (dB)');
        ylabel('Bit Error Rate (BER)');
        title(['BER Performance - SF = ' num2str(SF(j)) ', BW = ' num2str(BW(i)/1e6) ' MHz']);

        M = 2 .^ SF(j);
        t_TOA(i,j) = M ./ BW(i); % sec
        t_RS(i,j) = BW(i) ./ M;  % symbols per sec
    end

    
end