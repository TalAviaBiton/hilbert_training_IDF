close all; clc; clear;

%% data
delta_t = 0.1; % sec
sigma = 5; % meter per axis
sigma_a = 0.7; % meters / s ^ 2

%% section a: matrices
zeros4 = zeros(1, 4);
I = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

F = [1 delta_t 0 0; 0 1 0 0; 0 0 1 delta_t; 0 0 0 1]; 
H = [1 0 0 0; zeros4; 0 0 1 0; zeros4];
Q = sigma_a ^ 2 .* [delta_t ^ 4 / 4, delta_t ^ 3 / 2, delta_t ^ 3 / 2, delta_t ^ 2; zeros4; zeros4; zeros4]; 
R = [sigma 0 0 0; 0 0 0 0; 0 0 sigma 0; 0 0 0 0] .^ 2;

%% section b: course simulation

num_of_steps = 600;
steps_per_sec = num_of_steps * delta_t;

observation = zeros(num_of_steps, 4); % x = [Px, Vx, Py, Vy].';

a_x = sigma_a * randn(length(observation) - 1, 1);
a_y = sigma_a * randn(length(observation) - 1, 1);

observation(:, 2) = 12;
observation(:, 4) = 4;

for step = 2:num_of_steps

    observation(step, 2) = observation(step - 1, 2) + a_x(step - 1) * delta_t;
    observation(step, 4) = observation(step - 1, 4) + a_y(step - 1) * delta_t;

    observation(step, 1) = observation(step - 1, 1) + observation(step - 1, 2) * delta_t;
    observation(step, 3) = observation(step - 1, 3) + observation(step - 1, 4) * delta_t;

    % % section g: loosing every 3rd sample
    % if mod(step, 3) == 0
    %     observation(step, :) = [0 0 0 0];
    % end

end

figure;
plot(observation(:,1), observation(:,3), '-o', 'LineWidth', 2, 'MarkerSize', 6);
title('observed course');   
xlabel('X');  
ylabel('Y');
grid on; 

noise_x = randn(length(observation) - 1, 1);
noise_y = randn(length(observation) - 1, 1);
observation(2:end, 1) = observation(2:end, 1) + sigma .* noise_x;
observation(2:end, 3) = observation(2:end, 3) + sigma .* noise_y;

%% section c: kalman filter

Xk_k = observation(1, 1:4).';
Pk_k = [0.5 0 0 0; 0.5 0 0 0; 0.5 0 0 0; 0.5 0 0 0].';
Vk = det(R) * randn(size(observation, 2), 1);
x = zeros(num_of_steps, 4); % x = [Px, Vx, Py, Vy].'
x(1, :) = Xk_k;
p = zeros(num_of_steps, 4); 
p(1, :) = Pk_k(1, :);

prediction = zeros(num_of_steps, 4);
Y = zeros(num_of_steps, 4);
Z = zeros(num_of_steps, 4);
NIS = zeros(num_of_steps - 1, 4);
NEES = zeros(num_of_steps, 4);

for step = 2:num_of_steps

% prediction
Xk_k1 = F * Xk_k;
Pk_k1 = F * Pk_k * F.' + Q;
prediction(step, :) = Xk_k1;

% update
Zk = H * Xk_k1 + Vk;
Z(step, :) = Zk;

yk = Zk - H * Xk_k1;
Y(step, :) = yk;

Sk = (H * Pk_k1).' * H.' + R;

NIS(step-1, :) = yk.' * (Sk \ yk);
NIS(isnan(NIS)) = 0;

e = Xk_k - Xk_k1;
NEES(step, :) = e.' * (Pk_k1 \ e);
NEES(isnan(NEES)) = 0;

Kk = (Pk_k1.' * H.') ./ Sk;
Kk(isnan(Kk)) = 0;

Xk_k = Xk_k1 + Kk * yk;

Pk_k = (I - (Kk * H)) * [Pk_k(1, :); zeros4; zeros4; zeros4];

x(step, :) = Xk_k;
p(step, :) = Pk_k(1, :);

end

result = (x .* p + observation .* (1 - p));

figure;
plot(prediction(:,1), prediction(:,3), '-o', 'LineWidth', 0.2, 'MarkerSize', 6);
title('predicted course');   
xlabel('X');  
ylabel('Y');
grid on;  

figure;
plot(result(:,1), result(:,3), '-o', 'LineWidth', 0.2, 'MarkerSize', 6);
title('real course');   
xlabel('X');  
ylabel('Y');
grid on; 

figure;
plot(Y(:,1), Y(:,3), '-o', 'LineWidth', 2, 'MarkerSize', 6);
title('error');   
xlabel('X');  
ylabel('Y');
grid on;  

%% section d: RMSE

RMSE_observation = rmse(observation(:, 1:2:end), x(:, 1:2:end), "all");
RMSE_prediction = rmse(prediction(:, 1:2:end), x(:, 1:2:end), "all");
RMSE_velosity = rmse(prediction(:, 2:2:end), x(:, 2:2:end), "all");

%% section e: NIS & NEES

% NIS
figure;
plot(NIS, 'b-', 'LineWidth', 1.5, 'Color', "r");
hold on;
plot(x, 'b-', 'LineWidth', 1.5); 
title('Normalized Innovation Squared (NIS)');
xlabel('Time Step');
ylabel('NIS');
grid on;
avg_NIS = [mean(NIS(:,1), "all") mean(NIS(:,3), "all")];

% NEES
figure;
plot(NEES, 'm-', 'LineWidth', 1.5, 'Color', "r");
hold on;
plot(x, 'b-', 'LineWidth', 1.5); 
title('Normalized Estimation Error Squared (NEES)');
xlabel('Time Step');
ylabel('NEES');
grid on;
avg_NEES = [mean(NEES(:,1), "all") mean(NEES(:,3), "all")];


