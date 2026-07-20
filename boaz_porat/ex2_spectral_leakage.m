<<<<<<< HEAD
%% a
fc = 20e3;
fs = 1e6;
N = 100;
ts = 1/fs;
t = 0: ts:ts*(N-1);
x = sin( 2 * pi* t * fc);

%% b
f = -fs/2:fs/N:fs/2 - fs/N;
figure;
plot(t, x);
figure;
plot(f, fftshift(abs(fft(x))));

%% c
N = 120;
t = 0: ts:ts*(N-1);
x = sin( 2 * pi* t * fc);

f = -fs/2:fs/N:fs/2 - fs/N;
figure;
plot(t, x);
figure;
plot(f, fftshift(abs(fft(x))));
 % time resolution improved because of more samples wich leads to worse
=======
%% a
fc = 20e3;
fs = 1e6;
N = 100;
ts = 1/fs;
t = 0: ts:ts*(N-1);
x = sin( 2 * pi* t * fc);

%% b
f = -fs/2:fs/N:fs/2 - fs/N;
figure;
plot(t, x);
figure;
plot(f, fftshift(abs(fft(x))));

%% c
N = 120;
t = 0: ts:ts*(N-1);
x = sin( 2 * pi* t * fc);

f = -fs/2:fs/N:fs/2 - fs/N;
figure;
plot(t, x);
figure;
plot(f, fftshift(abs(fft(x))));
 % time resolution improved because of more samples wich leads to worse
>>>>>>> 8e4b46b2f4d6f53dcc1100746b1a1229c03dbe85
 % frequncy resolution