% =========================================================================
% Dynamic Tag Localization Simulation via Single Moving GPS Anchor (Tag A)
% Incorporating:
%   1. Markov-Chain Switchable Line-of-Sight / Non-Line-of-Sight Channel
%   2. Sage-Husa Adaptive Measurement Noise (R) Estimation
%   3. Huber M-Estimator Robust Residual Downweighting
% =========================================================================
clear; clc; close all;

%% 1. Simulation Setup & Configuration
dt = 0.01;                     % Sampling interval (seconds)
T_total = 15;                  % Total simulation duration (seconds)
t = 0:dt:T_total;
N = length(t);

%% 2. Ground Truth Motion Models (2D Trajectories)
% Tag A: Moving Leader (Anchor with GPS)
posA_true = [2*t; 5*sin(0.5*t)];               

% Tag B: Target Tag to be Estimated
posB_true = [1 + 1.5*t; -2 + 4*cos(0.3*t)];    

%% 3. Sensor Noise Parameters & Channel Simulation
sigma_gps  = 1.2;              % GPS noise standard deviation (meters)
sigma_los  = 0.10;             % UWB LOS noise standard deviation (meters)
mu_nlos    = 1.5;              % Mean exponential NLOS distance bias (meters)
sigma_nlos = 0.40;             % UWB NLOS noise standard deviation (meters)

% Generate Noisy Observations
posA_gps = posA_true + sigma_gps * randn(2, N);
true_dist = sqrt(sum((posA_true - posB_true).^2, 1));

% Synthesize UWB channel with Markov-chain driven NLOS multipath fading
[dist_uwb, nlos_states] = generate_uwb_nlos(true_dist, sigma_los, mu_nlos, sigma_nlos);

%% 4. Filter Setup (Extended Kalman Filter)
% State Vector: x = [x_pos; y_pos; x_vel; y_vel]
x_est = [0; 0; 0; 0];                          % Initial state guess
P = diag([10, 10, 2, 2]);                      % Initial error covariance
Q = diag([0.05, 0.05, 0.2, 0.2]);              % Process noise covariance
F = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];    % State transition matrix

% Robust & Adaptive Tuning Parameters
b = 0.96;                                      % Sage-Husa forgetting factor (0.95 - 0.99)
R_min = sigma_los^2;                           % Lower bound clamp for estimated R
R_est = R_min;                                 % Initial adaptive R estimate
huber_c = 2.0;                                 % Huber tuning threshold (standard deviations)

% Logging Allocation
posB_est    = zeros(2, N);
R_history   = zeros(1, N);
nis_history = zeros(1, N);

%% 5. Main Simulation Loop
for k = 1:N
    % --- Step 1: Time Prediction ---
    x_pred = F * x_est;
    P_pred = F * P * F' + Q;
    
    % --- Step 2: Measurement Prediction & Innovation ---
    pA = posA_gps(:, k);
    dx = x_pred(1) - pA(1);
    dy = x_pred(2) - pA(2);
    dist_pred = sqrt(dx^2 + dy^2);
    
    if dist_pred > 1e-4
        % Measurement Jacobian Matrix H = dh/dx
        H = [dx/dist_pred, dy/dist_pred, 0, 0];
        y = dist_uwb(k) - dist_pred;             % Innovation residual
        
        % --- Step 3: Sage-Husa Adaptive R Update ---
        d_k = (1 - b) / (1 - b^k);
        R_sample = (y^2) - (H * P_pred * H');
        R_est = (1 - d_k) * R_est + d_k * R_sample;
        R_est = max(R_est, R_min);               % Enforce floor bound
        
        % --- Step 4: Huber M-Estimator Weighting ---
        S_nominal = H * P_pred * H' + R_est;
        std_residual = abs(y) / sqrt(S_nominal);
        
        if std_residual <= huber_c
            w = 1.0;                             % Nominal weight
        else
            w = huber_c / std_residual;          % Downweight large residuals
        end
        
        % Inflate measurement variance robustly
        R_robust = R_est / w;
        S_robust = H * P_pred * H' + R_robust;
        
        % --- Step 5: Measurement Correction ---
        K = (P_pred * H') / S_robust;            % Kalman Gain
        x_est = x_pred + K * y;
        P = (eye(4) - K * H) * P_pred;
        
        % Log metrics
        nis_history(k) = (y^2) / S_robust;
    else
        x_est = x_pred;
        P = P_pred;
    end
    
    R_history(k) = R_est;
    posB_est(:, k) = x_est(1:2);
end

%% 6. Results & Visualization
rmse = sqrt(mean(sum((posB_true - posB_est).^2, 1)));
fprintf('===========================================\n');
fprintf('  Simulation Execution Complete\n');
fprintf('  Tag B Position Tracking RMSE: %.3f meters\n', rmse);
fprintf('===========================================\n');

figure('Name', 'Robust Adaptive UWB Localization', 'Position', [100, 100, 1100, 700]);

% Plot 1: 2D Spatial Trajectories
subplot(2, 2, [1, 3]);
plot(posA_true(1,:), posA_true(2,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Tag A Ground Truth'); hold on;
plot(posA_gps(1,:), posA_gps(2,:), 'c.', 'MarkerSize', 3, 'DisplayName', 'Tag A Noisy GPS');
plot(posB_true(1,:), posB_true(2,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Tag B Ground Truth');
plot(posB_est(1,:), posB_est(2,:), 'r:', 'LineWidth', 2, 'DisplayName', 'Tag B Robust EKF');
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('2D Target Tracking Trajectory');
legend('Location', 'best'); grid on; axis equal;

% Plot 2: UWB Ranging Distance
subplot(2, 2, 2);
plot(t, true_dist, 'g-', 'LineWidth', 1.5, 'DisplayName', 'True Distance'); hold on;
plot(t, dist_uwb, 'r.', 'MarkerSize', 4, 'DisplayName', 'Corrupted UWB Range');
xlabel('Time (s)'); ylabel('Distance (m)');
title('UWB Range Measurements (LOS vs NLOS Spikes)');
legend('Location', 'best'); grid on;

% Plot 3: Adaptive Covariance R Estimation
subplot(2, 2, 4);
plot(t, R_history, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Sage-Husa Estimated R'); hold on;
area(t, nlos_states * max(R_history)*0.8, 'FaceColor', [1 0.8 0.8], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'NLOS State Active');
xlabel('Time (s)'); ylabel('Covariance R (m^2)');
title('Real-Time Noise Covariance Estimation');
legend('Location', 'best'); grid on;

%% Local Function: Markov NLOS Channel Model
function [r_uwb, S_state] = generate_uwb_nlos(true_dist, sigma_los, mu_nlos, sigma_nlos)
    N_pts = length(true_dist);
    S_state = zeros(1, N_pts);
    r_uwb = zeros(1, N_pts);
    
    % Discrete-Time Markov Transition Matrix
    % [P(LOS->LOS), P(LOS->NLOS); P(NLOS->LOS), P(NLOS->NLOS)]
    P_trans = [0.98, 0.02; 
               0.10, 0.90];
    
    S_state(1) = 0; % Initial state: LOS
    for i = 1:N_pts
        if i > 1
            p_next = P_trans(S_state(i-1) + 1, :);
            S_state(i) = double(rand() > p_next(1));
        end
        
        if S_state(i) == 0  % Line-of-Sight (LOS)
            bias = 0;
            noise = sigma_los * randn();
        else              % Non-Line-of-Sight (NLOS)
            bias = exprnd(mu_nlos);
            noise = sigma_nlos * randn();
        end
        
        r_uwb(i) = true_dist(i) + bias + noise;
    end
end