%% Multi-Drone Navigation, EKF Sensor Fusion & Localization Error Analysis
% Author: AI Collaborator
% Description: Simulates Drone A (Leader) and Drone B (Follower) with EKF,
% ideal barometers, adaptive sensor reliability, and advanced error tracking.

clear; close all; clc;

%% 1. CONFIGURABLE SIMULATION PARAMETERS
config.coordSystem      = 'polar'; % Options: 'Cartesian' or 'Polar'
config.separationDist   = 30.0;         % Desired following distance (meters)
config.flightSpeed      = 1.0;         % Operational speed (m/s)
config.dt               = 0.1;         % Time step (seconds)
config.T_total          = 30.0;        % Total simulation time (seconds)

timeVec = 0:config.dt:config.T_total;
numSteps = length(timeVec);

%% 2. TRAJECTORY GENERATION FOR LEADER (DRONE A)
if strcmpi(config.coordSystem, 'Cartesian')
    omega = config.flightSpeed / 10.0; 
    true_xA = 10 * cos(omega * timeVec);
    true_yA = 10 * sin(omega * timeVec);
    true_zA = 2.0 + 0.1 * timeVec;     
else
    r = 10 + 0.05 * timeVec;
    theta = 0.5 * timeVec;
    true_xA = r .* cos(theta);
    true_yA = r .* sin(theta);
    true_zA = 2.0 + 0.05 * timeVec;
end

true_vxA = [diff(true_xA)/config.dt, 0];
true_vyA = [diff(true_yA)/config.dt, 0];

%% 3. TRAJECTORY GENERATION FOR FOLLOWER (DRONE B)
timeDelay = config.separationDist / config.flightSpeed;
delaySteps = round(timeDelay / config.dt);

true_xB = zeros(1, numSteps);
true_yB = zeros(1, numSteps);
true_zB = true_zA; 

for k = 1:numSteps
    refIndex = max(1, k - delaySteps);
    dx = true_xA(refIndex) - true_xA(k);
    dy = true_yA(refIndex) - true_yA(k);
    distNorm = sqrt(dx^2 + dy^2);
    if distNorm > 0
        true_xB(k) = true_xA(refIndex) - (config.separationDist * (dx / distNorm));
        true_yB(k) = true_yA(refIndex) - (config.separationDist * (dy / distNorm));
    else
        true_xB(k) = true_xA(k) - config.separationDist;
        true_yB(k) = true_yA(k);
    end
end

%% 4. SENSOR NOISE & INITIALIZATION
sigma_gps  = 1.5;   % GPS position noise (m)
sigma_uwb  = 0.5;   % UWB relative distance noise (m)
sigma_baro = 0.0;   % Barometer is 100% ideal / noise-free
sigma_dr_v = 0.2;   % Dead reckoning velocity noise (m/s)

est_A = zeros(6, numSteps);
est_B = zeros(6, numSteps);

P_A = eye(6) * 1.0;
P_B = eye(6) * 1.0;

reliability_A = ones(1, numSteps);
reliability_B = ones(1, numSteps);

%% 5. MAIN SIMULATION & EKF LOOP
for k = 2:numSteps
    dt = config.dt;
    
    %% --- DRONE A (LEADER) ---
    meas_gps_x  = true_xA(k) + randn * sigma_gps;
    meas_gps_y  = true_yA(k) + randn * sigma_gps;
    meas_baro_z = true_zA(k); 
    meas_dr_vx  = true_vxA(k) + randn * sigma_dr_v;
    meas_dr_vy  = true_vyA(k) + randn * sigma_dr_v;
    
    prev_xA = est_A(:, k-1);
    pred_x  = prev_xA(1) + prev_xA(4) * dt;
    pred_y  = prev_xA(2) + prev_xA(5) * dt;
    pred_z  = meas_baro_z; 
    x_pred_A = [pred_x; pred_y; pred_z; meas_dr_vx; meas_dr_vy; prev_xA(6)];
    
    baro_innovation = meas_baro_z - pred_z; 
    adaptive_gps_weight = 1.0 / (1.0 + abs(baro_innovation)); 
    reliability_A(k) = adaptive_gps_weight;
    
    H_A = eye(3, 6);
    R_A = diag([(sigma_gps / sqrt(adaptive_gps_weight))^2, ...
                (sigma_gps / sqrt(adaptive_gps_weight))^2, ...
                sigma_baro^2]);
    z_meas_A = [meas_gps_x; meas_gps_y; meas_baro_z];
    
    S_A = H_A * P_A * H_A' + R_A;
    K_A = P_A * H_A' / S_A;
    est_A(:, k) = x_pred_A + K_A * (z_meas_A - H_A * x_pred_A);
    P_A = (eye(6) - K_A * H_A) * P_A;

    %% --- DRONE B (FOLLOWER - GPS DENIED) ---
    true_rel_dist = sqrt((true_xA(k) - true_xB(k))^2 + (true_yA(k) - true_yB(k))^2);
    meas_uwb    = true_rel_dist + randn * sigma_uwb;
    meas_baro_zB = true_zB(k); 
    meas_dr_vxB = (true_xB(k) - true_xB(k-1))/dt + randn * sigma_dr_v;
    meas_dr_vyB = (true_yB(k) - true_yB(k-1))/dt + randn * sigma_dr_v;
    
    prev_xB = est_B(:, k-1);
    pred_xB_x = prev_xB(1) + prev_xB(4) * dt;
    pred_xB_y = prev_xB(2) + prev_xB(5) * dt;
    pred_xB_z = meas_baro_zB;
    x_pred_B = [pred_xB_x; pred_xB_y; pred_xB_z; meas_dr_vxB; meas_dr_vyB; prev_xB(6)];
    
    est_leader_pos = est_A(1:2, k);
    baro_innovation_B = meas_baro_zB - pred_xB_z;
    adaptive_uwb_weight = 1.0 / (1.0 + abs(baro_innovation_B));
    reliability_B(k) = adaptive_uwb_weight;
    
    H_B = eye(3, 6); 
    R_B = diag([(sigma_uwb / sqrt(adaptive_uwb_weight))^2, ...
                (sigma_uwb / sqrt(adaptive_uwb_weight))^2, ...
                sigma_baro^2]);
            
    z_meas_B = [est_leader_pos(1) - config.separationDist * cos(prev_xB(6)); ...
                est_leader_pos(2) - config.separationDist * sin(prev_xB(6)); ...
                meas_baro_zB];
            
    S_B = H_B * P_B * H_B' + R_B;
    K_B = P_B * H_B' / S_B;
    est_B(:, k) = x_pred_B + K_B * (z_meas_B - H_B * x_pred_B);
    P_B = (eye(6) - K_B * H_B) * P_B;
end

%% 6. VISUALIZATION & PLOTTING

% --- FIGURE 1: 3D Trajectory Tracking & Separation Performance ---
figure('Name', 'Flight Trajectories & Separation', 'Position', [100, 100, 1000, 600]);

subplot(1, 2, 1);
plot3(true_xA, true_yA, true_zA, 'b-', 'LineWidth', 2); hold on;
plot3(true_xB, true_yB, true_zB, 'r-', 'LineWidth', 2);
plot3(est_A(1,:), est_A(2,:), est_A(3,:), 'b--', 'LineWidth', 1.5);
plot3(est_B(1,:), est_B(2,:), est_B(3,:), 'r--', 'LineWidth', 1.5);
grid on; axis equal;
title(sprintf('3D Trajectory Tracking (%s)', config.coordSystem));
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Altitude Z (m)');
legend('True Leader', 'True Follower', 'EKF Leader', 'EKF Follower', 'Location', 'best');
view(50, 30);

subplot(1, 2, 2);
actual_separation = sqrt((true_xA - true_xB).^2 + (true_yA - true_yB).^2);
est_separation = sqrt((est_A(1,:) - est_B(1,:)).^2 + (est_A(2,:) - est_B(2,:)).^2);
plot(timeVec, actual_separation, 'k-', 'LineWidth', 1.5); hold on;
plot(timeVec, est_separation, 'm--', 'LineWidth', 1.5);
yline(config.separationDist, 'r:', 'Target Separation', 'LineWidth', 1.5);
grid on;
title('Inter-Drone Separation Distance');
xlabel('Time (s)'); ylabel('Distance (m)');
legend('True Separation', 'Estimated Separation', 'Target', 'Location', 'best');

% --- FIGURE 2: Drone B Localization Error & Sensor Metrics ---
figure('Name', 'Drone B Localization Error & Reliability', 'Position', [150, 150, 1100, 700]);

% Compute Drone B Position Errors
err_xB = est_B(1,:) - true_xB;
err_yB = est_B(2,:) - true_yB;
err_zB = est_B(3,:) - true_zB;
euclidean_err_B = sqrt(err_xB.^2 + err_yB.^2 + err_zB.^2);

subplot(2, 2, 1);
plot(timeVec, err_xB, 'r', 'LineWidth', 1.2); hold on;
plot(timeVec, err_yB, 'g', 'LineWidth', 1.2);
plot(timeVec, err_zB, 'b', 'LineWidth', 1.2);
grid on;
title('Drone B Error Components');
xlabel('Time (s)'); ylabel('Error (m)');
legend('X Error', 'Y Error', 'Z Error (Baro)', 'Location', 'best');

subplot(2, 2, 2);
plot(timeVec, euclidean_err_B, 'k-', 'LineWidth', 1.5);
grid on;
title('Drone B Total 3D Euclidean Localization Error');
xlabel('Time (s)'); ylabel('Total Error Magnitude (m)');

subplot(2, 2, [3, 4]);
plot(timeVec, reliability_A, 'b-', 'LineWidth', 1.5); hold on;
plot(timeVec, reliability_B, 'r-', 'LineWidth', 1.5);
grid on;
title('Adaptive Sensor Reliability Metric (Baro-Inferred Cross-Evaluation)');
xlabel('Time (s)'); ylabel('Reliability Weight (0 to 1)');
legend('Leader Reliability Factor', 'Follower Reliability Factor', 'Location', 'best');
ylim([-0.1, 1.2]);