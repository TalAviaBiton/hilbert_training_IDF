%% =========================================================================
% MULTI-DRONE NAVIGATION AND SENSOR FUSION SIMULATION (LEADER-FOLLOWER)
% =========================================================================
% Description: 
%   - Drone A (Leader): Sinusoidal forward flight path, GPS, Ideal Baro, DR.
%   - Drone B (Follower): GPS-denied, Ideal Baro, DR, and UWB relative ranging.
%   - Feature Added: Toggleable UWB calculation in Polar Coordinates (Radius 
%     from UWB range, angle estimation via EKF).
% =========================================================================

clear; close all; clc;

%% 1. CONFIGURABLE SIMULATION PARAMETERS
config.coord_system       = 'polar'; % Global display: 'cartesian' or 'polar'
config.uwb_representation = 'polar';     % UWB mode: 'polar' (radius calculated directly) or 'cartesian'
config.sep_distance       = 50.0;         % Following distance for Drone B (meters)
config.flight_speed       = 10.0;         % Operational speed (m/s)
config.dt                 = 0.001;         % Time step (seconds)
config.T_total            = 40.0;        % Total simulation time (seconds)

time_vec = 0:config.dt:config.T_total;
num_steps = length(time_vec);

%% 2. TRUE TRAJECTORY GENERATION (LEADER - DRONE A)
amplitude_y = 5.0;  % Sway amplitude (meters)
freq_y = 0.15;      % Sway frequency (rad/s)

true_A = struct();
true_A.x = zeros(1, num_steps);
true_A.y = zeros(1, num_steps);
true_A.z = zeros(1, num_steps);
true_A.v = zeros(1, num_steps);
true_A.psi = zeros(1, num_steps);

for k = 1:num_steps
    t = time_vec(k);
    true_A.x(k) = config.flight_speed * t;
    true_A.y(k) = amplitude_y * sin(freq_y * t);
    true_A.z(k) = 5.0; % Constant clean altitude profile
    
    if k == 1
        true_A.v(k) = config.flight_speed;
        true_A.psi(k) = 0;
    else
        dx = true_A.x(k) - true_A.x(k-1);
        dy = true_A.y(k) - true_A.y(k-1);
        true_A.v(k) = sqrt(dx^2 + dy^2) / config.dt;
        true_A.psi(k) = atan2(dy, dx);
    end
end

% Generate Follower (Drone B) Trajectory: lagging Drone A by separation distance
true_B = struct();
true_B.x = zeros(1, num_steps);
true_B.y = zeros(1, num_steps);
true_B.z = zeros(1, num_steps);
true_B.v = zeros(1, num_steps);
true_B.psi = zeros(1, num_steps);

for k = 1:num_steps
    lag_steps = round((config.sep_distance / config.flight_speed) / config.dt);
    ref_idx = max(1, k - lag_steps);
    
    true_B.x(k) = true_A.x(ref_idx) - config.sep_distance * cos(true_A.psi(ref_idx));
    true_B.y(k) = true_A.y(ref_idx) - config.sep_distance * sin(true_A.psi(ref_idx));
    true_B.z(k) = true_A.z(ref_idx);
    true_B.v(k) = true_A.v(ref_idx);
    true_B.psi(k) = true_A.psi(ref_idx);
end

%% 3. SENSOR NOISE & MEASUREMENT GENERATION
rng(42); % Seed for reproducibility

gps_noise_std = 1.5;     
dr_v_noise_std = 0.2;    
dr_psi_noise_std = 0.05; 
uwb_noise_std = 0.3;     

% Drone A Measurements
meas_A.gps_x = true_A.x + gps_noise_std * randn(1, num_steps);
meas_A.gps_y = true_A.y + gps_noise_std * randn(1, num_steps);
meas_A.baro_z = true_A.z; 
meas_A.dr_v = true_A.v + dr_v_noise_std * randn(1, num_steps);
meas_A.dr_psi = true_A.psi + dr_psi_noise_std * randn(1, num_steps);

% Drone B Measurements (GPS-Denied)
true_dist = sqrt((true_A.x - true_B.x).^2 + (true_A.y - true_B.y).^2 + (true_A.z - true_B.z).^2);
meas_B.uwb_dist = true_dist + uwb_noise_std * randn(1, num_steps);
meas_B.baro_z = true_B.z; 
meas_B.shared_alt_A = true_A.z; 
meas_B.dr_v = true_B.v + dr_v_noise_std * randn(1, num_steps);
meas_B.dr_psi = true_B.psi + dr_psi_noise_std * randn(1, num_steps);

%% 4. EXTENDED KALMAN FILTER (EKF) INITIALIZATION
est_A = zeros(5, num_steps);
est_A(:, 1) = [true_A.x(1); true_A.y(1); true_A.z(1); config.flight_speed; true_A.psi(1)];
P_A = eye(5) * 1.0;

est_B = zeros(5, num_steps);
est_B(:, 1) = [true_B.x(1); true_B.y(1); true_B.z(1); config.flight_speed; true_B.psi(1)];
P_B = eye(5) * 1.0;

% Storage for EKF Covariance (to track uncertainty bounds for Drone B)
P_history_B = zeros(5, num_steps);
P_history_B(:, 1) = diag(P_B);

Q = diag([0.01, 0.01, 0.01, 0.1, 0.05]);

%% 5. EKF ESTIMATION LOOP
for k = 2:num_steps
    %% --- DRONE A EKF UPDATE ---
    xA = est_A(1, k-1); yA = est_A(2, k-1); zA = est_A(3, k-1);
    vA = est_A(4, k-1); psiA = est_A(5, k-1);
    
    xA_pred = xA + vA * cos(psiA) * config.dt;
    yA_pred = yA + vA * sin(psiA) * config.dt;
    zA_pred = zA; 
    vA_pred = meas_A.dr_v(k);
    psiA_pred = meas_A.dr_psi(k);
    
    x_pred_A = [xA_pred; yA_pred; zA_pred; vA_pred; psiA_pred];
    
    F_A = [1, 0, 0, cos(psiA)*config.dt, -vA*sin(psiA)*config.dt;
           0, 1, 0, sin(psiA)*config.dt,  vA*cos(psiA)*config.dt;
           0, 0, 1, 0,            0;
           0, 0, 0, 1,            0;
           0, 0, 0, 0,            1];
       
    P_pred_A = F_A * P_A * F_A' + Q;
    
    z_meas_A = [meas_A.gps_x(k); meas_A.gps_y(k); meas_A.baro_z(k)];
    H_A = [1, 0, 0, 0, 0;
           0, 1, 0, 0, 0;
           0, 0, 1, 0, 0];
    R_A = diag([gps_noise_std^2, gps_noise_std^2, 1e-4]);
    
    y_res_A = z_meas_A - H_A * x_pred_A;
    S_A = H_A * P_pred_A * H_A' + R_A;
    K_A = P_pred_A * H_A' / S_A;
    
    est_A(:, k) = x_pred_A + K_A * y_res_A;
    P_A = (eye(5) - K_A * H_A) * P_pred_A;

    %% --- DRONE B EKF UPDATE (GPS-Denied) ---
    xB = est_B(1, k-1); yB = est_B(2, k-1); zB = est_B(3, k-1);
    vB = est_B(4, k-1); psiB = est_B(5, k-1);
    
    xB_pred = xB + vB * cos(psiB) * config.dt;
    yB_pred = yB + vB * sin(psiB) * config.dt;
    zB_pred = meas_B.shared_alt_A(k); 
    vB_pred = meas_B.dr_v(k);
    psiB_pred = meas_B.dr_psi(k);
    
    x_pred_B = [xB_pred; yB_pred; zB_pred; vB_pred; psiB_pred];
    
    F_B = [1, 0, 0, cos(psiB)*config.dt, -vB*sin(psiB)*config.dt;
           0, 1, 0, sin(psiB)*config.dt,  vB*cos(psiB)*config.dt;
           0, 0, 1, 0,            0;
           0, 0, 0, 1,            0;
           0, 0, 0, 0,            1];
       
    P_pred_B = F_B * P_B * F_B' + Q;
    
    % --- UWB MEASUREMENT PROCESSING: POLAR VS CARTESIAN OPTION ---
    dx_est = est_A(1, k) - xB_pred;
    dy_est = est_A(2, k) - yB_pred;
    dz_est = est_A(3, k) - zB_pred;
    
    if strcmpi(config.uwb_representation, 'polar')
        % Polar UWB formulation: Radius (range) is fully supplied by UWB,
        % while relative angle (bearing) and position are estimated by EKF.
        r_uwb_meas = meas_B.uwb_dist(k);
        r_est = sqrt(dx_est^2 + dy_est^2 + dz_est^2);
        
        z_meas_B = [r_uwb_meas; meas_B.shared_alt_A(k)];
        z_pred_B = [r_est; zB_pred];
        
        % Jacobian mapping polar radius and altitude to EKF state
        H_B = [-dx_est/r_est, -dy_est/r_est, -dz_est/r_est, 0, 0;
                0,              0,             1,             0, 0];
    else
        % Standard Cartesian distance formulation
        dist_est = sqrt(dx_est^2 + dy_est^2 + dz_est^2);
        
        z_meas_B = [meas_B.uwb_dist(k); meas_B.shared_alt_A(k)];
        z_pred_B = [dist_est; zB_pred];
        
        H_B = [-dx_est/dist_est, -dy_est/dist_est, -dz_est/dist_est, 0, 0;
                0,                0,                1,                  0, 0];
    end
    
    R_B = diag([uwb_noise_std^2, 1e-4]);
    y_res_B = z_meas_B - z_pred_B;
    
    S_B = H_B * P_pred_B * H_B' + R_B;
    K_B = P_pred_B * H_B' / S_B;
    
    est_B(:, k) = x_pred_B + K_B * y_res_B;
    P_B = (eye(5) - K_B * H_B) * P_pred_B;
    
    % Store covariance diagonal for Drone B uncertainty bounds
    P_history_B(:, k) = diag(P_B);
end

%% 6. VISUALIZATION SCRIPTS
% --- Figure 1: Flight Paths & Altitude Cross-Evaluation ---
figure('Name', 'Cartesian Trajectories & Sensor Fusion', 'Position', [100, 100, 1200, 550], 'Color', 'w');

subplot(1, 2, 1);
plot(true_A.x, true_A.y, 'b-', 'LineWidth', 2); hold on;
plot(est_A(1,:), est_A(2,:), 'b--', 'LineWidth', 1.5);
plot(true_B.x, true_B.y, 'r-', 'LineWidth', 2);
plot(est_B(1,:), est_B(2,:), 'r--', 'LineWidth', 1.5);
grid on; xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Leader (True)', 'Leader (EKF)', 'Follower (True)', 'Follower (EKF)', 'Location', 'best');
title(sprintf('Flight Path (UWB Mode: %s, Sep: %.1fm)', upper(config.uwb_representation), config.sep_distance));
axis equal;

subplot(1, 2, 2);
plot(time_vec, true_A.z, 'b-', 'LineWidth', 1.5); hold on;
plot(time_vec, est_A(3,:), 'b--', 'LineWidth', 1.2);
plot(time_vec, est_B(3,:), 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Altitude (m)');
legend('Leader True Alt', 'Leader EKF Alt', 'Follower EKF Alt (Shared Baro)', 'Location', 'best');
title('Altitude Sensor Fusion & Cross-Evaluation');

% --- Figure 2: Dedicated Drone B Localization Error & Uncertainty Bounds ---
figure('Name', 'Drone B Localization Error Analysis', 'Position', [150, 150, 1100, 700], 'Color', 'w');

err_xB = true_B.x - est_B(1,:);
err_yB = true_B.y - est_B(2,:);
err_zB = true_B.z - est_B(3,:);
err_total_B = sqrt(err_xB.^2 + err_yB.^2 + err_zB.^2);

sigma3_xB = 3 * sqrt(P_history_B(1,:));
sigma3_yB = 3 * sqrt(P_history_B(2,:));
sigma3_zB = 3 * sqrt(P_history_B(3,:));

subplot(3, 1, 1);
plot(time_vec, err_xB, 'r-', 'LineWidth', 1.5); hold on;
plot(time_vec,  sigma3_xB, 'k--', 'LineWidth', 1.2);
plot(time_vec, -sigma3_xB, 'k--', 'LineWidth', 1.2);
grid on; ylabel('X Error (m)');
legend('Actual X Error', '+/- 3\sigma EKF Bound', 'Location', 'best');
title(sprintf('Drone B Localization Errors & EKF Bounds (UWB Mode: %s)', upper(config.uwb_representation)));

subplot(3, 1, 2);
plot(time_vec, err_yB, 'r-', 'LineWidth', 1.5); hold on;
plot(time_vec,  sigma3_yB, 'k--', 'LineWidth', 1.2);
plot(time_vec, -sigma3_yB, 'k--', 'LineWidth', 1.2);
grid on; ylabel('Y Error (m)');
legend('Actual Y Error', '+/- 3\sigma EKF Bound', 'Location', 'best');

subplot(3, 1, 3);
plot(time_vec, err_total_B, 'b-', 'LineWidth', 2); hold on;
grid on; xlabel('Time (s)'); ylabel('3D Error (m)');
legend('Total 3D Position Error (Drone B)', 'Location', 'best');
title('Drone B Total 3D Position Localization Error');