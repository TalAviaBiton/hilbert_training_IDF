clc; clear; close all;

%% ex 1- freq resolution VS interpolation (zero padding)
clc; clear; close all;

f1=100;
f2=105;
Fs=1000;
phase=0;
Ts=1/Fs;
T=1;
t=0:Ts:T*(1-Ts);
f100=-Fs/2:Fs/100:Fs*(1/2-1/100);
f1000=-Fs/2:Fs/1000:Fs*(1/2-1/1000);
f1024=-Fs/2:Fs/1024:Fs*(1/2-1/1024);
signal=sin(2*pi*f1*t+phase)+sin(2*pi*f2*t+phase);
signal_100=signal(1:10:end);
fft_signal_100=fft(signal_100);
fft_signal_100_zero_padded=fft(signal_100,1024);
fft_signal_1000=fft(signal);

figure
plot(f100, abs(fft_signal_100));
xlabel('freq (Hz)');
ylabel('Amplitude');
title('fft signal 100');

figure
plot(f1024, abs(fft_signal_100_zero_padded));
xlabel('freq (Hz)');
ylabel('Amplitude');
title('fft signal 100 zero padded');

figure
plot(f1000, abs(fft_signal_1000));
xlabel('freq (Hz)');
ylabel('Amplitude');
title('fft signal 1000');

%% 2 windowing
clc; clear; close all;

fc=100.5;
phase=0;
Fs=100;
Ts=1/Fs;
L=10;
T=1;
N=100;

t=0:Ts:T*(1-Ts);
f=-Fs/2:Fs/N:Fs*(1/2-1/N);
signal=sin(2*pi*fc*t+phase);

rect_sig=pulse_signal([ones(1),zeros(1)], signal);
hamming_sig=pulse_signal(hamming(L), signal);
blackman_sig=pulse_signal(blackman(L), signal);

display(t,f,N,signal,'no window');
display(t,f,N,rect_sig,'rect');
display(t,f,N,hamming_sig,'hamming');
display(t,f,N,blackman_sig,'blackman');

%% 3 spectrogram
clc; clear; close all;

T=5;
Fs=200000; %time resolution
Ts=1/Fs;
t=0:Ts:T*(1-Ts);
f0=1000;
t1=50;
f1=50000;
signal=chirp(t,f0,t1,f1);
% figure
% spectrogram(signal);

L=512; %freq resolution
win=hamming(L);
overlap=L*0.25; %time resolution
nfft=1024; %freq resolution
[spec,freq,time,power] = spectrogram(signal,win,overlap,nfft,Fs);

figure;
imagesc(time, freq, 10*log10(power)); 
axis xy;                        
colormap(jet);                  
colorbar;                       
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Signal Spectrogram');

%% 4 Downsampling and LPF
clc; clear; close all;

Fc=400;
Fs=1000;
phase=0;
T=1;
jumps=4;
N=Fs/jumps;
Ts=1/Fs;
t=0:Ts:T*(1-Ts);
f=-Fs/2:Fs/N:Fs*(1/2-1/N);
f_og=-Fs/2:Fs/(N*jumps):Fs*(1/2-1/(N*jumps));

signal=sin(2*pi*Fc*t+phase);
filterd_sig=lowpass(signal,Fs/jumps, Fs);
filterd_sig_ds=filterd_sig(1:jumps:end);
samples_sig1=signal(1:jumps:end);
samples_sig2=decimate(signal,jumps);

display(t,f_og,N*jumps,signal,'og')
display(t,f_og,N*jumps,filterd_sig,'filterd')
display(t,f,N,filterd_sig_ds,'filterd ds')
display(t,f,N,samples_sig1,'every 4th sample')
display(t,f,N,samples_sig2,'built in decimation function')

%% 5 correlation
clc; clear; close all;

% n=100;%ceil(1/(rand())); 
% signal=randi([0,1],1,n);
fc=4;
phase=0;
Fs=100;
Ts=1/Fs;
L=10;
T=1;
N=100;

t=0:Ts:T*(1-Ts);
f=-Fs/2:Fs/N:Fs*(1/2-1/N);
signal=sin(2*pi*fc*t+phase);

shift=randi([0,1],1,ceil(1/(rand()))); 
shifted_signal=[shift,signal];

noise=randi([0,1],1,length(shifted_signal));
noised_signal=shifted_signal+noise;

padded_signal=[zeros(1,length(noised_signal)-length(signal)),signal];

 display((0:length(signal)),f,N,signal,'og');
% display((0:length(shift)),f,N,shift,'shift');
% display((0:length(shifted_signal)),f,N,shifted_signal,'shifted');
% display((0:length(noise)),f,N,noise,'noise');
% display((0:length(noised_signal)),f,N,noised_signal,'noised');

[corelation,l]=xcorr(noised_signal,padded_signal);
[peak,indx]=max(corelation);
% figure
% plot(corelation);

decoded_shift_len=indx-length(signal);
deshifted_signal=noised_signal(decoded_shift_len:end);
decoded_signal=deshifted_signal.^2;
thrashold=1;
decoded_signal(decoded_signal>thrashold)=1;
decoded_signal(decoded_signal<thrashold)=0;

display((0:length(decoded_signal)),f,N,decoded_signal,'decoded');


%% functions


function pulsed_signal=pulse_signal(window, signal)
    
    sig_len=length(signal);
    win_len=length(window);
    reps=ceil(sig_len/win_len);
    rep_win=repmat(window,1,reps);
    pulse=rep_win(1:sig_len);
    pulsed_signal=pulse.*signal;
    
end

function display(t,f,N,signal,graph_title)
    
    jumps=length(t)/length(signal);

    figure
    plot(t(1:jumps:end), signal);
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title(graph_title );%+ ': time');

    % fft_signal=fft(signal,N);
    % 
    % figure
    % plot(f, abs(fft_signal));
    % xlabel('freq (Hz)');
    % ylabel('Amplitude');
    % title(graph_title );%+ ': freq');

end

function filterd_sig = lowpass_filter(signal,f)
fft_signal=fft(signal);
l=length(fft_signal)/2;
fft_signal(l+f+1:end)=0;
fft_signal(1:l-f-1)=0;
filterd_sig=ifft(fft_signal);
end
