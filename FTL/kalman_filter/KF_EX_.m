clc; clear;
%% data
delta_t = 0.1; % sec
sigma = 5; % meter per axis (measurement noise std)
sigma_a = 0.7; % meters / s ^ 2 (process noise std)

%% section a: matrices
zeros2 = zeros(2, 2);
I = eye(4);
F = [1 delta_t 0 0; 
     0 1      0 0; 
     0 0      1 delta_t; 
     0 0      0 1]; 

% We measure only Position X and Position Y -> H is 2x4
H = [1 0 0 0; 
     0 0 1 0];

% 1D process noise block for [p, v]
q_block = sigma_a^2 * [delta_t^4/4, delta_t^3/2; 
                       delta_t^3/2, delta_t^2];
% 2D decoupled Q matrix (4x4)
Q = [q_block, zeros2; 
     zeros2,  q_block]; 

% Measurement covariance R for [z_px, z_py] -> 2x2
R = (sigma^2) * eye(2);

%% section b: course simulation
num_of_steps = 600;
observation = zeros(num_of_steps, 4); % x = [Px, Vx, Py, Vy].'
a_x = sigma_a * randn(num_of_steps - 1, 1);
a_y = sigma_a * randn(num_of_steps - 1, 1);

observation(1, 2) = 12; % initial velocity X
observation(1, 4) = 4;  % initial velocity Y

for step = 2:num_of_steps
    observation(step, 2) = observation(step - 1, 2) + a_x(step - 1) * delta_t;
    observation(step, 4) = observation(step - 1, 4) + a_y(step - 1) * delta_t;
    observation(step, 1) = observation(step - 1, 1) + observation(step - 1, 2) * delta_t;
    observation(step, 3) = observation(step - 1, 3) + observation(step - 1, 4) * delta_t;
end

% Add measurement noise to positions
noise_x = randn(num_of_steps, 1);
noise_y = randn(num_of_steps, 1);
observed_positions = zeros(num_of_steps, 2);
observed_positions(:, 1) = observation(:, 1) + sigma .* noise_x;
observed_positions(:, 2) = observation(:, 3) + sigma .* noise_y;

%% section c: kalman filter
Xk_k = [observation(1, 1); observation(1, 2); observation(1, 3); observation(1, 4)];
Pk_k = eye(4) * 10; % Initial error covariance

x = zeros(num_of_steps, 4); 
x(1, :) = Xk_k';
prediction = zeros(num_of_steps, 4);
NIS = zeros(num_of_steps - 1, 1);
NEES = zeros(num_of_steps, 1);

for step = 2:num_of_steps
    % 1. Prediction step
    Xk_k1 = F * Xk_k;
    Pk_k1 = F * Pk_k * F' + Q;
    prediction(step, :) = Xk_k1';
    
    % 2. Update step
    Zk = observed_positions(step, :)'; % True measurement at step
    yk = Zk - H * Xk_k1;
    
    Sk = H * Pk_k1 * H' + R;
    
    % NIS (Normalized Innovation Squared)
    NIS(step-1) = yk' * (Sk \ yk);
    
    % Kalman Gain
    Kk = Pk_k1 * H' / Sk;
    
    % State update
    Xk_k = Xk_k1 + Kk * yk;
    
    % Covariance update (Joseph form)
    Pk_k = (I - Kk * H) * Pk_k1 * (I - Kk * H)' + Kk * R * Kk';
    
    % NEES (Normalized Estimation Error Squared)
    e = Xk_k - observation(step, :)';
    NEES(step) = e' * (Pk_k \ e);
    
    x(step, :) = Xk_k';
end

%% section d: RMSE
RMSE_position = rmse(observation(:, [1,3]), x(:, [1,3]), "all");
RMSE_velocity = rmse(observation(:, [2,4]), x(:, [2,4]), "all");

%% section e: NIS & NEES plots
figure;
plot(NIS, 'b-', 'LineWidth', 1.5);
title('Normalized Innovation Squared (NIS)');
xlabel('Time Step');
ylabel('NIS');
grid on;

figure;
plot(NEES, 'r-', 'LineWidth', 1.5);
title('Normalized Estimation Error Squared (NEES)');
xlabel('Time Step');
ylabel('NEES');
grid on;

figure;
plot(observation(:,1), observation(:,3), 'k--', 'LineWidth', 1.5);
hold on;
plot(x(:,1), x(:,3), 'r-', 'LineWidth', 1.5);
legend('True Trajectory', 'KF Estimate');
title('Real vs Estimated Course');
xlabel('X'); ylabel('Y');
grid on;