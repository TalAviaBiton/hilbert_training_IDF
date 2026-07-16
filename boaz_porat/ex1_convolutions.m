%% PART A
N = 30;
M = 50;

signal1 = symetryc_rect(N, 10);
signal2 = symetryc_rect(M, 20);

stem(conv(signal1,signal2,"full"));

function rect = symetryc_rect(len, num_of_samples)

pad = zeros(1,(len-num_of_samples)/2);
rect = [pad ones(1, num_of_samples) pad];

end

%% PART B
fs = 100; %Hz
fc = 20; %Hz
N = 1000;
SNR = 6; %dB
ts = 1/fs; %sec

t = 0:ts:ts * (N-1);
f = -fs/2:fs/N:fs/2 - fs/N;

sig = sin(2*pi*fc*t);
y = awgn(sig, SNR, "measured");
figure;
plot(t, (y));

% section 1
B = 40;%Hz
if B > fs/2
    out = 'error'
end
A = 2*B/fs;
n = (-N/2:N/2-1);
x = A *sinc(A*n);

% section 2
figure;
plot(t, abs(x));
%i expected a rect but i got almost rect
figure;
plot(f, fftshift(abs(fft(x))));

%section 3
%i expect to get a mess in the BW
stem(conv(x,y,"full"));
