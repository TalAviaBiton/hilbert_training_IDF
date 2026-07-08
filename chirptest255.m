% הגדרות (באותם ערכים מהשלב הקודם)
fs = 1000;           % העלינו את תדר הדגימה ל-1000 כדי לקבל רזולוציה טובה יותר בספקטוגרמה
t = 0:1/fs:0.5;      % זמן ארוך יותר להמחשה ברורה
f0 = 10;             % תדר התחלה
f1 = 400;            % תדר סיום
sig = chirp(t, f0, 0.5, f1);

% יצירת הספקטוגרמה
figure;
window = hamming(64);       % חלון ניתוח
noverlap = 32;              % חפיפה בין חלונות
nfft = 128;                 % נקודות ב-FFT

spectrogram(sig, window, noverlap, nfft, fs, 'yaxis');
title('Chirp Signal Spectrogram');

% נורמליזציה לטווח של 0-255 (עבור 8 ביטים)
sig_norm = round((sig - min(sig)) / (max(sig) - min(sig)) * 255);

% המרה לביטים (כל דגימה ל-8 ביטים)
bits = dec2bin(sig_norm, 8);
bit_stream = reshape(bits', 1, []);

% בדיקת אורך
fprintf('Total bits: %d\n', length(bit_stream));
disp('First 32 bits:');
disp(bit_stream(1:32));