close all; clc; clear;

x = (1:100);
N = 3;
polyphased_changed_fs_x = Polyphase_filter(x, N);
stem(polyphased_changed_fs_x);

function polyphased_changed_fs_x = Polyphase_filter(x, N)
    if mod(N, 1) == 0 || N < 1
        polyphased_changed_fs_x = change_fs(x, N);
    else
        [phases, ~] = rat(N);
        
        x_phases = reshape(x, [], phases);
        l = length(x);
        x_phases_changed = zeros(ceil(l * N / phases), phases);
        
        for i = 1:phases
            tmp = change_fs(x_phases(:, i), N).';
            x_phases_changed(:, i) = [tmp ; zeros(size(x_phases_changed, 1) - length(tmp), 1)];
        end
        polyphased_changed_fs_x = x_phases_changed(:).';
    end
end

function x_change_fs = change_fs(x, N)
    [L, M] = rat(N); 
    interpulated_x = InterpulatePP(x, L);
    x_change_fs = DecimatePP(interpulated_x, M);
end

function decimaed_x = DecimatePP(x, M)
    if M <= 1
        decimaed_x = x;
        return;
    end
    Wn = 1 / M;
    n = 40;
    filter_h = fir1(n, Wn); 
    
    v = conv(x, filter_h, 'same'); 
    decimaed_x = v(1:M:end);
end

function interpulated_x = InterpulatePP(x, L)
    if L <= 1 
        interpulated_x = x;
        return;
    end
    up_sampled_x = upsample(x, L);
    Wn = 1 / L;
    n = 40; 
    filter_h = L * fir1(n, Wn);
    
    interpulated_x = conv(up_sampled_x, filter_h, 'same');
end