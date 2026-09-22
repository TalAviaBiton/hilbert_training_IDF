close all; clc; clear;

x = (1:100);
N = 0.8; % עכשיו יעבוד לכל מספר שבר או שלם

polyphased_changed_fs_x = Polyphase_filter(x, N);

% נצייר את התוצאות כדי לראות את ההבדל
figure;
subplot(2,1,1);
stem(x);
title('Original Signal');

subplot(2,1,2);
stem(polyphased_changed_fs_x);
title(['Resampled Signal (N = ' num2str(N) ')']);

function y = Polyphase_filter(x, N)
    [L, M] = rat(N);
    
    if L == 1 && M == 1
        y = x;
        return;
    end
    
    x_up = upsample(x, L);
    
    Wn = min(1/L, 1/M);
    
    n = 60; 
    
    filter_h = L * fir1(n, Wn);
    
    v = conv(x_up, filter_h, 'same');
    
    y = v(1:M:end);
end