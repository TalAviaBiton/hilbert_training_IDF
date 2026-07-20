spec = create_spectrogram(signal, Fs, dt, RBW, overlap, rect);

function spec = create_spectrogram(signal, Fs, dt, RBW, overlap, window)
% dt - requierd time resolution
% RBW - requierd frequncy resolution
% overlap - precentage of overlap between two windows
% window - anonimic function that returns the window
buf_sig= buffer(signal, RBW, overlap);
win = window(length(signal), )
f_start = 1e6;
f_stop = 2e6;
t_start = 0;
t_stop = 1e3;
Fc = 1.5e6;

%check rbw and dt corespond to eachother
Ts = 1 / Fs;
t = t_start:1 / dt:t_stop - 1 / dt;
f = f_start:1 / RBW:f_stop - 1 / RBW;
imagesc(t, f, buf_sig);
set(buf_sig, 'ydir', 'normal')

end

function rect = symetryc_rect(len, num_of_samples)

pad = zeros(1,(len-num_of_samples)/2);
rect = [pad ones(1, num_of_samples) pad];

end