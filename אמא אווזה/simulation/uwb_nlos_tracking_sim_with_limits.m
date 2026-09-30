% =========================================================================
% 3D Dynamic Tag Localization Simulation via Kinematically Clamped Anchor
% Features: 3D Kinematics, Tag A Kinematic Bounds (Distance, Velocity, Acc),
%           Selectable Noise, Sage-Husa Adaptive R, Huber Robust EKF
% =========================================================================
clear; clc;% close all;

%% 1. Simulation Setup & Configuration
dt = 0.01;                     % Sampling interval (seconds)
T_total = 15;                  % Total simulation duration (seconds)
t = 0:dt:T_total;
N = length(t);

%% 2. Ground Truth 3D Motion Models [X; Y; Z]
posA_true = [2*t; 5*sin(0.5*t); 10 + 0.5*t];               
posB_true = [1 + 1.5*t; -2 + 4*cos(0.3*t); 5 + 2*sin(0.2*t)];    

%% 3. Sensor Noise & Tag A Kinematic Bounds
noise_type  = 'bursty';        % Choose: 'gaussian' or 'bursty'
p_burst     = 0.05;            % Burst occurrence probability
sigma_burst = 8.0;             % Burst noise std dev (meters)

sigma_gps   = 1.5;             % GPS noise standard deviation (meters)
sigma_los   = 0.10;            % UWB LOS noise standard deviation (meters)
mu_nlos     = 1.5;             % Mean exponential NLOS distance bias (meters)
sigma_nlos  = 0.40;            % UWB NLOS noise standard deviation (meters)

% --- TAG A KINEMATIC LIMITS STRUCTURE ---
lim.dist_min_axis = [0.001, 0.001, 0.0005];  % Min step distance per axis (m)
lim.dist_max_axis = [0.080, 0.080, 0.0300];  % Max step distance per axis (m)
lim.vel_min_axis  = [0.050, 0.050, 0.0200];  % Min velocity per axis (m/s)
lim.vel_max_axis  = [3.000, 3.000, 1.5000];  % Max velocity per axis (m/s)
lim.acc_min_axis  = [0.000, 0.000, 0.0000];  % Min acceleration per axis (m/s^2)
lim.acc_max_axis  = [2.000, 2.000, 1.0000];  % Max acceleration per axis (m/s^2)

lim.dist_min_total = 0.002;                  % Min 3D total step distance (m)
lim.dist_max_total = 0.100;                  % Max 3D total step distance (m)
lim.vel_min_total  = 0.100;                  % Min 3D total velocity (m/s)
lim.vel_max_total  = 4.000;                  % Max 3D total velocity (m/s)
lim.acc_min_total  = 0.000;                  % Min 3D total acceleration (m/s^2)
lim.acc_max_total  = 2.500;                  % Max 3D total acceleration (m/s^2)

% Generate Raw Noisy Observations
posA_gps_raw = posA_true + sigma_gps * randn(3, N);
posA_clamped = zeros(3, N);
velA_clamped = zeros(3, N);

% Apply Kinematic Clamping Pre-Filter to Tag A GPS Data
posA_clamped(:, 1) = posA_gps_raw(:, 1);
velA_clamped(:, 1) = [0; 0; 0];
for k = 2:N
    [posA_clamped(:, k), velA_clamped(:, k)] = clamp_kinematics(...
        posA_gps_raw(:, k), posA_clamped(:, k-1), velA_clamped(:, k-1), dt, lim);
end

% Synthesize UWB Channel based on Clamped Tag A Position
true_dist = sqrt(sum((posA_true - posB_true).^2, 1));
[dist_uwb, nlos_states, ~] = generate_uwb_channel_3d(...
    true_dist, sigma_los, mu_nlos, sigma_nlos, noise_type, p_burst, sigma_burst);

%% 4. Filter Setup (3D EKF for Tag B)
x_est = [0; 0; 0; 0; 0; 0];                    
P = diag([10, 10, 10, 2, 2, 2]);               
Q = diag([0.05, 0.05, 0.05, 0.2, 0.2, 0.2]);   
F = [eye(3), dt*eye(3); zeros(3), eye(3)];     

b = 0.96;                                      
R_min = sigma_los^2;                           
R_est = R_min;                                 
huber_c = 2.0;                                 

posB_est  = zeros(3, N);
R_history = zeros(1, N);

%% 5. Main Simulation Loop
for k = 1:N
    % Step 1: Time Prediction
    x_pred = F * x_est;
    P_pred = F * P * F' + Q;
    
    % Step 2: 3D Innovation using Kinematically Clamped Tag A Position
    pA = posA_clamped(:, k);
    dx = x_pred(1) - pA(1);
    dy = x_pred(2) - pA(2);
    dz = x_pred(3) - pA(3);
    dist_pred = sqrt(dx^2 + dy^2 + dz^2);
    
    if dist_pred > 1e-4
        H = [dx/dist_pred, dy/dist_pred, dz/dist_pred, 0, 0, 0];
        y = dist_uwb(k) - dist_pred;             
        
        % Step 3: Sage-Husa Adaptive R
        d_k = (1 - b) / (1 - b^k);
        R_sample = (y^2) - (H * P_pred * H');
        R_est = max((1 - d_k) * R_est + d_k * R_sample, R_min);
        
        % Step 4: Huber M-Estimator Weighting
        S_nominal = H * P_pred * H' + R_est;
        std_residual = abs(y) / sqrt(S_nominal);
        w = iff(std_residual <= huber_c, 1.0, huber_c / std_residual);
        
        S_robust = H * P_pred * H' + (R_est / w);
        
        % Step 5: Measurement Correction
        K = (P_pred * H') / S_robust;            
        x_est = x_pred + K * y;
        P = (eye(6) - K * H) * P_pred;
    else
        x_est = x_pred;
        P = P_pred;
    end
    
    R_history(k) = R_est;
    posB_est(:, k) = x_est(1:3);
end

%% 6. Visualization & Metrics
rmse_3d = sqrt(mean(sum((posB_true - posB_est).^2, 1)));
fprintf('===========================================\n');
fprintf('  Tag B 3D Tracking RMSE: %.3f meters\n', rmse_3d);
fprintf('===========================================\n');

figure('Name', '3D Tracking with Kinematic Clamping', 'Position', [100, 100, 1200, 700]);

subplot(1, 2, 1);
plot3(posA_true(1,:), posA_true(2,:), posA_true(3,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Tag A True'); hold on;
plot3(posA_gps_raw(1,:), posA_gps_raw(2,:), posA_gps_raw(3,:), 'c.', 'MarkerSize', 2, 'DisplayName', 'Tag A Raw GPS');
plot3(posA_clamped(1,:), posA_clamped(2,:), posA_clamped(3,:), 'm-', 'LineWidth', 1.5, 'DisplayName', 'Tag A Clamped GPS');
plot3(posB_true(1,:), posB_true(2,:), posB_true(3,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Tag B True');
plot3(posB_est(1,:), posB_est(2,:), posB_est(3,:), 'r:', 'LineWidth', 2, 'DisplayName', 'Tag B EKF');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('3D Trajectories with Kinematic Clamping');
legend('Location', 'best'); grid on; view(3);

subplot(1, 2, 2);
plot(t, sqrt(sum((posA_gps_raw - posA_true).^2, 1)), 'c-', 'DisplayName', 'Raw GPS Error'); hold on;
plot(t, sqrt(sum((posA_clamped - posA_true).^2, 1)), 'm-', 'LineWidth', 1.5, 'DisplayName', 'Clamped GPS Error');
xlabel('Time (s)'); ylabel('Position Error (m)');
title('Tag A Position Error Reduction via Kinematic Limits');
legend('Location', 'best'); grid on;

%% --- HELPER FUNCTIONS ---
function [pos_out, vel_out] = clamp_kinematics(pos_raw, pos_prev, vel_prev, dt, lim)
    d_raw = pos_raw - pos_prev;
    v_raw = d_raw / dt;
    a_raw = (v_raw - vel_prev) / dt;
    
    % 1. Acceleration Clamping (Axis & Total)
    a_clamped = clamp_vector(a_raw, lim.acc_min_axis, lim.acc_max_axis, lim.acc_min_total, lim.acc_max_total);
    
    % 2. Velocity Clamping (Axis & Total)
    v_target = vel_prev + a_clamped * dt;
    v_clamped = clamp_vector(v_target, lim.vel_min_axis, lim.vel_max_axis, lim.vel_min_total, lim.vel_max_total);
    
    % 3. Displacement Clamping (Axis & Total)
    d_target = v_clamped * dt;
    d_clamped = clamp_vector(d_target, lim.dist_min_axis, lim.dist_max_axis, lim.dist_min_total, lim.dist_max_total);
    
    pos_out = pos_prev + d_clamped;
    vel_out = d_clamped / dt;
end

function vec_out = clamp_vector(vec_in, min_axis, max_axis, min_total, max_total)
    vec_out = vec_in;
    
    % Per-axis clamping
    for i = 1:3
        mag = abs(vec_out(i));
        sgn = sign(vec_out(i));
        if sgn == 0, sgn = 1; end
        mag_clamped = min(max(mag, min_axis(i)), max_axis(i));
        vec_out(i) = sgn * mag_clamped;
    end
    
    % Total 3D Norm Clamping
    norm_val = norm(vec_out);
    if norm_val > 1e-9
        norm_clamped = min(max(norm_val, min_total), max_total);
        vec_out = vec_out * (norm_clamped / norm_val);
    end
end

function [r_uwb, S_state, burst_log] = generate_uwb_channel_3d(true_dist, sigma_los, mu_nlos, sigma_nlos, noise_type, p_burst, sigma_burst)
    N_pts = length(true_dist);
    S_state = zeros(1, N_pts); burst_log = zeros(1, N_pts); r_uwb = zeros(1, N_pts);
    P_trans = [0.98, 0.02; 0.10, 0.90]; S_state(1) = 0;
    for i = 1:N_pts
        if i > 1, S_state(i) = double(rand() > P_trans(S_state(i-1)+1, 1)); end
        bias = iff(S_state(i) == 0, 0, exprnd(mu_nlos));
        noise = iff(S_state(i) == 0, sigma_los * randn(), sigma_nlos * randn());
        if strcmpi(noise_type, 'bursty') && (rand() < p_burst)
            burst_log(i) = sigma_burst * randn();
        end
        r_uwb(i) = true_dist(i) + bias + noise + burst_log(i);
    end
end

function val = iff(cond, t_val, f_val)
    if cond, val = t_val; else, val = f_val; end
end