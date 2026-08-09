
x = (1:100);
M = 3;
fc = ;
fs = ;

function polyphased_changed_fs_x = Polyphase_filter(x, N, phases)

phases = mindivisor(N, phases);

x_phases = reshape(x, [], phases);

for i = 1:phases
    change_fs(x_phases(i), N);
end

polyphased_changed_fs_x = x_phases(:).';

end

function x_change_fs = change_fs (x, N)

l = length(x);
M = mindivisor(l, 0);
L = N / M;

interpulated_x = InterpulatePP (x, L);
x_change_fs = DecimatePP (interpulated_x, M);

end

function decimaed_x = DecimatePP (x, M)

m1 = 1;
m2 = length(x);
Len = ; % dependent on n and m2 maybe more
extended_x = Extend (x, Len, m1, m2);

rad = ;
sample = ;
Wn = (1 / M) * (pi * rad / sample);
n = ;
filter = fir1(n, Wn);

v = conv(extended_x, filter);

y = dis_Extend (v, m1, m2);

decimaed_x = y(1:M:end);

decimaed_x = decimaed_x(1: floor(l *M));

end

function interpulated_x = InterpulatePP (x, L)

up_sampled_x = upsample(x, L);

m1 = 1;
m2 = length(up_sampled_x);
Len = ; % dependent on n and m2 maybe more
extended_x = Extend (up_sampled_x, Len, m1, m2);

rad = ;
sample = ;
Wn = (1 / L) * (pi * rad / sample);
n = ;
filter = fir1(n, Wn);

v = conv(extended_x, filter);

interpulated_x = dis_Extend (v, m1, m2);

end

function extended_x = Extend (x, L, m1, m2)

    for i = m1 - L:m1
        left_extantion = 2 * x(m1) - x(2 * m1 - n);
    end
    
    for i = m2:m2 + L
        right_extantion = 2 * x(m2) - x(2 * m2 - n);
    end
    
    extended_x = [left_extantion x right_extantion];

end

function x = dis_Extend (extended_x, m1, m2)

    x = extended_x(m1:m2);

end

function m = mindivisor(N, a)
    divs = N ./ (1:N);
    for i = 1:N
        if mod(divs(i), 1) ~= 0
            divs(i) = inf(1,1);
        end
    end
    [~, idx] = min(abs(divs - a));
    m= divs(idx);
end
