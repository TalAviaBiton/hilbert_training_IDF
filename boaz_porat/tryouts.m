close all;

fs = 50;             % קצב דגימה בהרץ
ts = 1/fs
t = 0:ts:1;            % וקטור זמן למשך שנייה אחת
f = 5;                 % תדר אות הסינוס בהרץ
amplitude = 1;         % משרעת (אמפליטודה)

y = amplitude * sin(2 * pi * f * t); % חישוב האות

plot(t, y);            % הצגת הגרף
xlabel('זמן (שניות)');
ylabel('משרעת');
title('אות סינוס');
grid on;
