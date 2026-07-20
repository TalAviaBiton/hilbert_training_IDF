
%%% Example %%%
rf_freq = 470e6;
sf = 11;
bw = 1.25e6;
fs = 1e6;
phy = LoRaPHY(rf_freq, sf, bw, fs);
phy.has_header = 1; % explicit header mode

symbols = [2541,1153,673,2397,1189,3509,41,3089,3237,3917,2729,2765,1417,2833,1389,801,3197,345,961,745,3101,297,1893,469]';
[data, checksum] = phy.decode(symbols);
disp(data); % CODE: 09 90 40 01 02 03 04 05 06 07 08 09 BA 2E
disp(checksum);

