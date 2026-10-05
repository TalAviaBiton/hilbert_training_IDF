% =========================================================================
% 3D Dynamic Tag Localization: UWB + Dual Barometer Fusion (With Error Graph)
% Features: Tag B moves with identical velocity profile to Tag A,
%           Constant 100m separation
% A has barometer & gps and B has only the barometer. B tracks itself based
% on UWB with A.
% =========================================================================
clear; clc; close all;

%% 1. Configuration
dt = 0.01;                     % Sampling interval (seconds)
T_total = 100;                  % Total simulation duration (seconds)
t = 0:dt:T_total;
N = length(t);

%% 2. Ground Truth 3D Motion Models

distance_betwwen_A_and_B = 30; % meters

% Tag A: Moving Leader (Anchor)
xA = 2*t;
yA = 3*sin(0.3*t);
zA = 0.5*t;
posA_true = [xA; yA; zA];               

% Tag B: Follower Tag - (Identical velocities)
posB_true = posA_true - [distance_betwwen_A_and_B; 0; 0];    

%% 3. Sensor Noise & Barometer Configuration
noise_type  = 'bursty';        % Choose: 'gaussian' or 'bursty'
p_burst     = 0.05;            % Burst occurrence probability (5%)
sigma_burst = 8.0;             % High-amplitude burst noise std dev (meters)

sigma_gps   = 1.2;             % GPS horizontal noise standard deviation (meters)
sigma_baro  = 0.15;            % Barometer height noise standard deviation (meters)
sigma_los   = 0.10;            % UWB LOS noise standard deviation (meters)
mu_nlos     = 1.5;             % Mean exponential NLOS distance bias (meters)
sigma_nlos  = 0.40;            % UWB NLOS noise standard deviation (meters)

% Generate Noisy Sensor Observations
posA_gps = posA_true + sigma_gps * randn(3, N);
baroA_meas = posA_true(3, :) + sigma_baro * randn(1, N);
baroB_meas = posB_true(3, :) + sigma_baro * randn(1, N);

% Fuse Tag A's Barometer into its Z coordinate
posA_gps(3, :) = baroA_meas; 

% True UWB Distance between Tag A and Tag B (Constant 100m)
true_dist = sqrt(sum((posA_true - posB_true).^2, 1));
[dist_uwb, nlos_states, burst_log] = generate_uwb_channel_3d(...
    true_dist, sigma_los, mu_nlos, sigma_nlos, noise_type, p_burst, sigma_burst);

%% 4. EKF setup
% --- FILTER TYPE SELECTION ---
filter_type = 'polar';         % Choose: 'cartesian' or 'polar'

% Adaptive & Robust Tuning Parameters
b = 0.96;                      % Sage-Husa forgetting factor
R_min_uwb = sigma_los^2;       % Lower bound clamp for UWB R
R_uwb_est = R_min_uwb;         
huber_c = 2.0;                 % Huber tuning threshold

F = [eye(3), dt*eye(3); zeros(3), eye(3)]; % 6x6 State transition matrix

if strcmpi(filter_type, 'cartesian')
    % Absolute Cartesian State: x = [x; y; z; vx; vy; vz]
    x_est = [posB_true(:, 1) + [2; -1; 0.5]; 0; 0; 0];  
    P = diag([10, 10, 5, 2, 2, 1]);
    Q = diag([0.05, 0.05, 0.01, 0.2, 0.2, 0.05]);
else
    % Relative Polar State: x = [r; theta; phi; r_dot; theta_dot; phi_dot]
    pA_0 = posA_gps(:, 1);
    pB_0 = posB_true(:, 1);
    rel_0 = pB_0 - pA_0; 
    r_0     = norm(rel_0);     % Exactly 100 meters
    theta_0 = atan2(rel_0(2), rel_0(1));
    phi_0   = atan2(rel_0(3), sqrt(rel_0(1)^2 + rel_0(2)^2));
    
    x_est = [r_0; theta_0; phi_0; 0; 0; 0];
    P = diag([5, 0.5, 0.5, 1, 0.1, 0.1]);
    Q = diag([0.02, 0.001, 0.001, 0.1, 0.01, 0.01]);
end

% Logging Allocation
posB_est    = zeros(3, N);
R_history   = zeros(1, N);

%% 5. Main Simulation Loop (Multi-Sensor Fusion: UWB + Barometer)
for k = 1:N
    % --- Step 1: Time Prediction ---
    x_pred = F * x_est;
    P_pred = F * P * F' + Q;
    
    pA = posA_gps(:, k);       
    zA = pA(3);                
    zB_baro = baroB_meas(k);   
    z_uwb = dist_uwb(k);       
    
    z_meas = [z_uwb; zB_baro];
    
    % --- Step 2: Measurement Prediction & Jacobian Setup ---
    if strcmpi(filter_type, 'cartesian')
        dx = x_pred(1) - pA(1);
        dy = x_pred(2) - pA(2);
        dz = x_pred(3) - pA(3);
        dist_pred = sqrt(dx^2 + dy^2 + dz^2);
        
        z_pred = [dist_pred; x_pred(3)]; 
        
        if dist_pred > 1e-4 % zero devising prevention
            H = [dx/dist_pred, dy/dist_pred, dz/dist_pred, 0, 0, 0;
                 0,            0,            1,            0, 0, 0];
        else
            H = zeros(2, 6);
        end
    else
        r_p = x_pred(1);
        th_p = x_pred(2);
        ph_p = x_pred(3);
        
        dist_pred = r_p;
        height_pred = zA + r_p * sin(ph_p);
        
        z_pred = [dist_pred; height_pred];
        
        if r_p > 1e-4
            H = [1,  0,  0,               0, 0, 0;
                 sin(ph_p), 0, r_p*cos(ph_p), 0, 0, 0];
        else
            H = zeros(2, 6);
        end
    end
    
    y = z_meas - z_pred; 
    
    if dist_pred > 1e-4
        % --- Step 3: Sage-Husa Adaptive R for UWB ---
        d_k = (1 - b) / (1 - b^k);
        R_uwb_sample = (y(1)^2) - (H(1,:) * P_pred * H(1,:)');
        R_uwb_est = (1 - d_k) * R_uwb_est + d_k * R_uwb_sample;
        R_uwb_est = max(R_uwb_est, R_min_uwb);       
        
        % --- Step 4: Huber M-Estimator Weighting on UWB ---
        S_uwb_nom = H(1,:) * P_pred * H(1,:)' + R_uwb_est;
        std_residual = abs(y(1)) / sqrt(S_uwb_nom);
        
        if std_residual <= huber_c
            w = 1.0; 
        else
            w = huber_c / std_residual;          
        end
        
        R_uwb_robust = R_uwb_est / w;
        
        R_mat = diag([R_uwb_robust, sigma_baro^2]);
        S_robust = H * P_pred * H' + R_mat;
        
        % --- Step 5: Measurement Correction ---
        K = (P_pred * H') / S_robust;            
        x_est = x_pred + K * y;
        P = (eye(6) - K * H) * P_pred;
    else
        x_est = x_pred;
        P = P_pred;
    end
    
    % --- Step 6: Map Estimated State back to Cartesian Position for Output ---
    if strcmpi(filter_type, 'cartesian')
        posB_est(:, k) = x_est(1:3);
    else
        r_e     = x_est(1);
        theta_e = x_est(2);
        phi_e   = x_est(3);
        
        dx_est = r_e * cos(phi_e) * cos(theta_e);
        dy_est = r_e * cos(phi_e) * sin(theta_e);
        dz_est = r_e * sin(phi_e);
        
        posB_est(:, k) = pA + [dx_est; dy_est; dz_est];
    end
    
    R_history(k) = R_uwb_est;
end

%% 6. Results & Visualization
rmse_3d = sqrt(mean(sum((posB_true - posB_est).^2, 1)));
fprintf('===========================================\n');
fprintf('  Selected Filter Mode: %s EKF\n', upper(filter_type));
fprintf('  Selected Noise Mode : %s\n', upper(noise_type));
fprintf('  Tag B 3D Position Tracking RMSE: %.3f meters\n', rmse_3d);
fprintf('===========================================\n');

% Calculate 3D Euclidean error over time
error_3d_vector = sqrt(sum((posB_true - posB_est).^2, 1));
mean_error = mean(error_3d_vector);

% Plot 1: 3D Spatial Trajectories (Top-Left)
figure
plot3(posA_true(1,:), posA_true(2,:), posA_true(3,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Tag A (Leader)'); hold on;
plot3(posB_true(1,:), posB_true(2,:), posB_true(3,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Tag B (Truth)');
plot3(posB_est(1,:), posB_est(2,:), posB_est(3,:), 'r:', 'LineWidth', 2, 'DisplayName', 'Tag B (Estimated)');
xlabel('X Position (m)'); ylabel('Y Position (m)'); zlabel('Z Position (m)');
title('3D Trajectory Tracking');
legend('Location', 'best'); grid on; view(3);

% Plot 2: Tag B Localization Error over Time (Top-Right)
figure
plot(t, error_3d_vector, 'r-', 'LineWidth', 1.5, 'DisplayName', '3D Error Magnitude'); hold on;
yline(mean_error, 'k--', sprintf('Mean Error: %.2f m', mean_error), 'LineWidth', 1.5, 'DisplayName', 'Mean Error');
xlabel('Time (s)'); ylabel('Error (m)');
title('Tag B Independent Localization Error');
legend('Location', 'best'); grid on;

%% Local Functions
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
        
        if S_state(i) == 0  
            bias = 0;
            noise = sigma_los * randn();
        else              
            bias = exprnd(mu_nlos);
            noise = sigma_nlos * randn();
        end
        
        if strcmpi(noise_type, 'bursty')
            if rand() < p_burst
                burst_log(i) = sigma_burst * randn(); 
            end
        end
        
        r_uwb(i) = true_dist(i) + bias + noise + burst_log(i);
    end
end