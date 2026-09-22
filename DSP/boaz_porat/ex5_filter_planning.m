%% definitions

Ap = 1; % גליות בתחום המעבר בdB
As = 40; % ניחות בתחום החסימה dB-
fp = 475; % תדר מעבר Hz
fs = 500; % תדר חסימה Hz
Fs = 100e3; % תדר דגימה מקורי Hz
M = 100; % decimation factor

%% part a - conversion of A to delta(linear) and f to theta(radians)

deltaP = (10^(Ap/20) - 1) / (10^(Ap/20) + 1);
deltaS = 10^(-As/20);

thetaP = 2 * pi * fp;
thetaS = 2 * pi * fs;

%% part b - plan a LPF, plot its frequency response and find the filter order (N).

[n, fo, ao, w] = firpmord([fp, fs], [1, 0], [deltaP, deltaS], Fs);
LPF = firpm(n, fo, ao, w);
l = length(LPF);
f = -Fs/2:Fs/l:Fs/2 - Fs/l;

figure
plot(f, abs(LPF));
xlabel('freq (Hz)');
ylabel('Amplitude');
title(' frequency response of the LPF');

N = n;

MPS = (N + 1) / 60;

%% part c - change Ap As and fs - fp (transition band width) and diagnose the affects on N, which has the largest impact?

%% part d - divide the system to two step decimation: 20 then 5. (same As & Ap). find the filters order with firpmord. how many multiplication operations are needed per second?