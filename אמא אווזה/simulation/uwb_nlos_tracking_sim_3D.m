% =========================================================================
% 3D Dynamic Tag Localization Simulation via Single Moving GPS Anchor (Tag A)
% Features: 3D Kinematics, Selectable Gaussian / Bursty Channel Noise,
%           Sage-Husa Adaptive R, Huber Robust EKF
% =========================================================================
clear; clc;% close all;

%% 1. Simulation Setup & Configuration
dt = 0.01;                     % Sampling interval (seconds)
T_total = 15;                  % Total simulation duration (seconds)
t = 0:dt:T_total;
N = length(t);

%% 2. Ground Truth 3D Motion Models [X; Y; Z]
% Tag A: Moving Leader (Anchor with 3D GPS)
posA_true = [2*t; 5*sin(0.5*t); 10 + 0.5*t];               

% Tag B: Target Tag to be Estimated
posB_true = [1 + 1.5*t; -2 + 4*cos(0.3*t); 5 + 2*sin(0.2*t)];    

%% 3. Sensor Noise & Selectable Channel Model Parameters
% --- NOISE TYPE SELECTION ---
noise_type  = 'bursty';        % Choose: 'gaussian' or 'bursty'

% Bursty Noise Parameters (Bernoulli-Gaussian Model)
p_burst     = 0.05;            % Burst occurrence probability (5% of measurements)
sigma_burst = 8.0;             % High-amplitude burst noise std dev (meters)

% Standard Sensor Noise Parameters
sigma_gps   = 1.2;             % GPS noise standard deviation (meters)
sigma_los   = 0.10;            % UWB LOS noise standard deviation (meters)
mu_nlos     = 1.5;             % Mean exponential NLOS distance bias (meters)
sigma_nlos  = 0.40;            % UWB NLOS noise standard deviation (meters)

% Generate Noisy Observations
posA_gps = posA_true + sigma_gps * randn(3, N);
true_dist = sqrt(sum((posA_true - posB_true).^2, 1));

% Synthesize UWB Channel (LOS, NLOS, and optional Bursty Noise)
[dist_uwb, nlos_states, burst_log] = generate_uwb_channel_3d(...
    true_dist, sigma_los, mu_nlos, sigma_nlos, noise_type, p_burst, sigma_burst);

%% 4. Filter Setup (3D Extended Kalman Filter)
x_est = [0; 0; 0; 0; 0; 0];                    % Initial state guess [x; y; z; vx; vy; vz]
P = diag([10, 10, 10, 2, 2, 2]);               % Initial error covariance (6x6)
Q = diag([0.05, 0.05, 0.05, 0.2, 0.2, 0.2]);   % Process noise covariance (6x6)
F = [eye(3), dt*eye(3); zeros(3), eye(3)];     % 6x6 State transition matrix

% Adaptive & Robust Tuning Parameters
b = 0.96;                                      % Sage-Husa forgetting factor
R_min = sigma_los^2;                           % Lower bound clamp for estimated R
R_est = R_min;                                 % Initial adaptive R estimate
huber_c = 2.0;                                 % Huber tuning threshold

% Logging Allocation
posB_est    = zeros(3, N);
R_history   = zeros(1, N);
nis_history = zeros(1, N);

%% 5. Main Simulation Loop
for k = 1:N
    % --- Step 1: Time Prediction ---
    x_pred = F * x_est;
    P_pred = F * P * F' + Q;
    
    % --- Step 2: 3D Measurement Prediction & Innovation ---
    pA = posA_gps(:, k);
    dx = x_pred(1) - pA(1);
    dy = x_pred(2) - pA(2);
    dz = x_pred(3) - pA(3);
    dist_pred = sqrt(dx^2 + dy^2 + dz^2);
    
    if dist_pred > 1e-4
        H = [dx/dist_pred, dy/dist_pred, dz/dist_pred, 0, 0, 0];
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
            w = huber_c / std_residual;          % Downweight large bursts/outliers
        end
        
        % Inflate measurement variance robustly
        R_robust = R_est / w;
        S_robust = H * P_pred * H' + R_robust;
        
        % --- Step 5: Measurement Correction ---
        K = (P_pred * H') / S_robust;            % Kalman Gain (6x1)
        x_est = x_pred + K * y;
        P = (eye(6) - K * H) * P_pred;
        
        nis_history(k) = (y^2) / S_robust;
    else
        x_est = x_pred;
        P = P_pred;
    end
    
    R_history(k) = R_est;
    posB_est(:, k) = x_est(1:3);
end

%% 6. Results & Visualization
rmse_3d = sqrt(mean(sum((posB_true - posB_est).^2, 1)));
fprintf('===========================================\n');
fprintf('  Selected Noise Mode: %s\n', upper(noise_type));
fprintf('  Tag B 3D Position Tracking RMSE: %.3f meters\n', rmse_3d);
fprintf('===========================================\n');

figure('Name', sprintf('3D EKF Tracking - [%s Noise]', upper(noise_type)), 'Position', [100, 100, 1200, 750]);

% Plot 1: 3D Spatial Trajectories
subplot(2, 2, [1, 3]);
plot3(posA_true(1,:), posA_true(2,:), posA_true(3,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Tag A Ground Truth'); hold on;
plot3(posA_gps(1,:), posA_gps(2,:), posA_gps(3,:), 'c.', 'MarkerSize', 3, 'DisplayName', 'Tag A Noisy GPS');
plot3(posB_true(1,:), posB_true(2,:), posB_true(3,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Tag B Ground Truth');
plot3(posB_est(1,:), posB_est(2,:), posB_est(3,:), 'r:', 'LineWidth', 2, 'DisplayName', 'Tag B Robust EKF');
xlabel('X Position (m)'); ylabel('Y Position (m)'); zlabel('Z Position (m)');
title(sprintf('3D Trajectory Tracking (%s Mode)', upper(noise_type)));
legend('Location', 'best'); grid on; view(3);

% Plot 2: UWB Range Measurements & Bursts
subplot(2, 2, 2);
plot(t, true_dist, 'g-', 'LineWidth', 1.5, 'DisplayName', 'True Distance'); hold on;
plot(t, dist_uwb, 'r.', 'MarkerSize', 4, 'DisplayName', 'Corrupted UWB Range');
if strcmpi(noise_type, 'bursty')
    burst_idx = find(burst_log ~= 0);
    plot(t(burst_idx), dist_uwb(burst_idx), 'ko', 'MarkerSize', 6, 'DisplayName', 'Bursts Detected');
end
xlabel('Time (s)'); ylabel('Distance (m)');
title('UWB Range Measurements');
legend('Location', 'best'); grid on;

% Plot 3: Adaptive Covariance R Estimation
subplot(2, 2, 4);
plot(t, R_history, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Sage-Husa Estimated R'); hold on;
area(t, nlos_states * max(R_history)*0.8, 'FaceColor', [1 0.8 0.8], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'NLOS State Active');
xlabel('Time (s)'); ylabel('Covariance R (m^2)');
title('Adaptive Measurement Noise R');
legend('Location', 'best'); grid on;

%% Local Function: 3D Channel Model with Selectable Bursty Noise
function [r_uwb, S_state, burst_log] = generate_uwb_channel_3d(...
    true_dist, sigma_los, mu_nlos, sigma_nlos, noise_type, p_burst, sigma_burst)

    N_pts = length(true_dist);
    S_state   = zeros(1, N_pts);
    burst_log = zeros(1, N_pts);
    r_uwb     = zeros(1, N_pts);
    
    P_trans = [0.98, 0.02; 
               0.10, 0.90];
    
    S_state(1) = 0;
    for i = 1:N_pts
        if i > 1
            p_next = P_trans(S_state(i-1) + 1, :);
            S_state(i) = double(rand() > p_next(1));
        end
        
        % Nominal Noise & NLOS Bias
        if S_state(i) == 0  % Line-of-Sight (LOS)
            bias = 0;
            noise = sigma_los * randn();
        else              % Non-Line-of-Sight (NLOS)
            bias = exprnd(mu_nlos);
            noise = sigma_nlos * randn();
        end
        
        % Optional Bursty Noise Injection (Bernoulli-Gaussian)
        if strcmpi(noise_type, 'bursty')
            if rand() < p_burst
                burst_log(i) = sigma_burst * randn(); % High-magnitude impulse
            end
        end
        
        r_uwb(i) = true_dist(i) + bias + noise + burst_log(i);
    end
end