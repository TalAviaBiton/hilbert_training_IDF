% close all;
clear; clc;

%% experiments
f1=1e6;
f2=0.25e6;
Fs=3e6;
phase=0;
Ts=1/Fs;
T=1;
t=0:Ts:T*(1-Ts);
sig=chirp(t,f1,100, f2);
bw=3e6;
sf=7;
specp(sig, Fs, bw, sf);



figure('Name', 'Spectrogram of Received LoRa Signal');
spectrogram(sig, 256, 250, 256, Fs, 'yaxis');
title('Spectrogram of Last Modulated Signal');

%% function spec

function specp(sig, fs, bw, sf)
x = 0:1/sf:2^sf/bw;
y = -bw/2:bw/2;
p = fs/bw;
win_length = 2^(sf-2);
N = round(p*2^sf);
[s] = spectrogram(sig,win_length,round(win_length*0.8),N);
valid_data_len = round(2^sf/2*1.5);
b = abs(s(valid_data_len:-1:1,:));
c = abs(s(end:-1:end-valid_data_len+1,:));
d = [b;c];
figure;
imagesc(x,y,d);
colormap summer;
title('Spectrogram');
xlabel('Time');
ylabel('Frequency');
set(gcf,'unit','normalized','position',[0.05,0.2,0.9,0.1]);
end