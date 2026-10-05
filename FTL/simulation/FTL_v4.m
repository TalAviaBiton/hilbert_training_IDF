%% =========================================================================
% MULTI-DRONE NAVIGATION & SENSOR FUSION SIMULATION (WITH MEAN ERROR LINES)
% =========================================================================
% Description: Drone A follows a sinusoidal trajectory (GPS + Dead Reckoning).
% At every delta step, Drone B calculates its distance (radius) and direction 
% (angle) from Drone A, using polar relative data for dynamic control via UWB 
% and an EKF. Includes mean error reference lines on error graphs.
% =========================================================================

clear; clc; close all;

%% 1. Configurable Simulation Parameters
config.coord_system        = 'Polar';     % State representation in EKF
config.separation_distance = 30;          % Target relative radius from Drone A (meters)
config.flight_speed        = 5;           % Operational flight speed (m/s)
config.dt                  = 0.1;         % Time step delta (seconds)
config.total_time          = 30;          % Total simulation time (seconds)

time_vector = 0:config.dt:config.total_time;
num_steps   = length(time_vector);

%% 2. Drone A (Leader) True Trajectory Generation (Sinusoidal Path)
x_A_true = config.flight_speed * time_vector;
y_A_true = 30 * sin(0.15 * time_vector);

vx_A_true = gradient(x_A_true, config.dt);
vy_A_true = gradient(y_A_true, config.dt);

r_A_true     = sqrt(x_A_true.^2 + y_A_true.^2);
theta_A_true = atan2(y_A_true, x_A_true);

%% 3. Sensor Noise Specifications
sigma_gps_r     = 1.5;   % GPS position error radius (meters)
sigma_gps_theta = 0.02;  % GPS position error angle (radians)

sigma_uwb_dist  = 0.5;   % UWB relative distance noise (meters)
sigma_uwb_angle = 0.01;  % UWB relative angle noise (radians)

%% 4. Extended Kalman Filter (EKF) Initialization
% State Vector: x = [r; theta; v_r; v_theta]

% Drone A EKF Initialization
x_est_A = [r_A_true(1); theta_A_true(1); vx_A_true(1); vy_A_true(1)];
P_A     = diag([1, 0.1, 1, 0.1]);
Q_A     = diag([0.01, 0.001, 0.1, 0.01]); 
R_A     = diag([sigma_gps_r^2, sigma_gps_theta^2]); 

% Pre-allocate Drone B True & Estimated Trajectories
x_B_true = zeros(1, num_steps);
y_B_true = zeros(1, num_steps);

% Initial position: Drone B starts 50m behind Drone A
x_B_true(1) = x_A_true(1) - config.separation_distance;
y_B_true(1) = y_A_true(1);

r_B_init     = sqrt(x_B_true(1)^2 + y_B_true(1)^2);
theta_B_init = atan2(y_B_true(1), x_B_true(1));

x_est_B = [r_B_init; theta_B_init; 0; 0];
P_B     = diag([1, 0.1, 1, 0.1]);
Q_B     = diag([0.02, 0.002, 0.2, 0.02]);
R_B     = diag([sigma_uwb_dist^2, sigma_uwb_angle^2]);

est_traj_A = zeros(4, num_steps);
est_traj_B = zeros(4, num_steps);

%% 5. Simulation Loop (Delta-Step Polar Control & EKF)
for k = 1:num_steps
    dt = config.dt;

    % =====================================================
    % DRONE B DYNAMIC PATH CONTROL (Calculated at every delta)
    % =====================================================
    if k > 1
        % 1. Calculate distance (radius) and direction (angle) from Drone A at this delta
        dx_rel = x_A_true(k-1) - x_B_true(k-1);
        dy_rel = y_A_true(k-1) - y_B_true(k-1);
        
        r_from_A     = sqrt(dx_rel^2 + dy_rel^2);       % Distance (Radius) from A
        angle_from_A = atan2(dy_rel, dx_rel);           % Direction (Angle) from A
        
        % 2. Determine next step direction and speed based on polar relative data
        radius_error = r_from_A - config.separation_distance;
        
        % Leader heading estimation
        if k > 2
            lead_heading = atan2(y_A_true(k-1) - y_A_true(k-2), x_A_true(k-1) - x_A_true(k-2));
        else
            lead_heading = 0;
        end
        
        % Proportional control to maintain separation distance and track leader heading
        speed_cmd = config.flight_speed + 1.2 * radius_error;
        speed_cmd = max(0, min(speed_cmd, 2 * config.flight_speed));
        
        % Direction vector towards formation target point behind Drone A
        target_x = x_A_true(k-1) - config.separation_distance * cos(lead_heading);
        target_y = y_A_true(k-1) - config.separation_distance * sin(lead_heading);
        
        dir_x = target_x - x_B_true(k-1);
        dir_y = target_y - y_B_true(k-1);
        dist_to_target = sqrt(dir_x^2 + dir_y^2);
        
        if dist_to_target > 0
            vx_B_true_k = speed_cmd * (dir_x / dist_to_target);
            vy_B_true_k = speed_cmd * (dir_y / dist_to_target);
        else
            vx_B_true_k = 0; vy_B_true_k = 0;
        end
        
        % Update Drone B true position for current delta
        x_B_true(k) = x_B_true(k-1) + vx_B_true_k * dt;
        y_B_true(k) = y_B_true(k-1) + vy_B_true_k * dt;
    end

    % =====================================================
    % DRONE A EKF (Leader: Dead Reckoning + GPS)
    % =====================================================
    r_pred_A     = x_est_A(1) + x_est_A(3) * dt;
    theta_pred_A = x_est_A(2) + x_est_A(4) * dt;
    vr_pred_A    = x_est_A(3);
    vtheta_pred_A= x_est_A(4);
    
    x_pred_A = [r_pred_A; theta_pred_A; vr_pred_A; vtheta_pred_A];
    F_A = [1, 0, dt, 0; 0, 1, 0, dt; 0, 0, 1, 0; 0, 0, 0, 1];
    P_pred_A = F_A * P_A * F_A' + Q_A;

    z_gps_r     = r_A_true(k) + sigma_gps_r * randn();
    z_gps_theta = theta_A_true(k) + sigma_gps_theta * randn();
    z_A         = [z_gps_r; z_gps_theta];

    H_A = [1, 0, 0, 0; 0, 1, 0, 0];
    y_A = z_A - H_A * x_pred_A;
    S_A = H_A * P_pred_A * H_A' + R_A;
    K_A = P_pred_A * H_A' / S_A;
    
    x_est_A = x_pred_A + K_A * y_A;
    P_A     = (eye(4) - K_A * H_A) * P_pred_A;
    est_traj_A(:, k) = x_est_A;

    % =====================================================
    % DRONE B EKF (Follower: Dead Reckoning + UWB Relative)
    % =====================================================
    r_pred_B     = x_est_B(1) + x_est_B(3) * dt;
    theta_pred_B = x_est_B(2) + x_est_B(4) * dt;
    vr_pred_B    = x_est_B(3);
    vtheta_pred_B= x_est_B(4);
    
    x_pred_B = [r_pred_B; theta_pred_B; vr_pred_B; vtheta_pred_B];
    F_B      = F_A;
    P_pred_B = F_B * P_B * F_B' + Q_B;

    % True relative UWB measurements (Distance & Angle relative to Drone A)
    r_B_true_k     = sqrt(x_B_true(k)^2 + y_B_true(k)^2);
    theta_B_true_k = atan2(y_B_true(k), x_B_true(k));

    xA_k = r_A_true(k) * cos(theta_A_true(k));
    yA_k = r_A_true(k) * sin(theta_A_true(k));
    xB_k = r_B_true_k * cos(theta_B_true_k);
    yB_k = r_B_true_k * sin(theta_B_true_k);

    true_rel_dist  = sqrt((xB_k - xA_k)^2 + (yB_k - yA_k)^2);
    true_rel_angle = atan2(yB_k - yA_k, xB_k - xA_k);

    z_uwb_dist  = true_rel_dist + sigma_uwb_dist * randn();
    z_uwb_angle = true_rel_angle + sigma_uwb_angle * randn();
    z_B         = [z_uwb_dist; z_uwb_angle];

    % Measurement model based on estimated states
    r_b = x_pred_B(1); th_b = x_pred_B(2);
    r_a = x_est_A(1);  th_a = x_est_A(2); 

    xB_est = r_b * cos(th_b); yB_est = r_b * sin(th_b);
    xA_est = r_a * cos(th_a); yA_est = r_a * sin(th_a);

    dx_rel        = xB_est - xA_est;
    dy_rel        = yB_est - yA_est;
    pred_rel_dist = sqrt(dx_rel^2 + dy_rel^2);
    pred_rel_angle= atan2(dy_rel, dx_rel);

    h_B = [pred_rel_dist; pred_rel_angle];

    ddist_drb  = (dx_rel * cos(th_b) + dy_rel * sin(th_b)) / pred_rel_dist;
    ddist_dthb = (dx_rel * (-r_b * sin(th_b)) + dy_rel * (r_b * cos(th_b))) / pred_rel_dist;
    dangle_drb = (-dy_rel * cos(th_b) + dx_rel * sin(th_b)) / (pred_rel_dist^2);
    dangle_dthb= (-dy_rel * (-r_b * sin(th_b)) + dx_rel * (r_b * cos(th_b))) / (pred_rel_dist^2);

    H_B = [ddist_drb,   ddist_dthb,   0, 0;
           dangle_drb,  dangle_dthb,  0, 0];

    y_B = z_B - h_B;
    y_B(2) = atan2(sin(y_B(2)), cos(y_B(2)));

    S_B = H_B * P_pred_B * H_B' + R_B;
    K_B = P_pred_B * H_B' / S_B;

    x_est_B = x_pred_B + K_B * y_B;
    P_B     = (eye(4) - K_B * H_B) * P_pred_B;
    est_traj_B(:, k) = x_est_B;
end

%% 6. Post-Processing & Error Calculation
x_B_est_cart = est_traj_B(1, :) .* cos(est_traj_B(2, :));
y_B_est_cart = est_traj_B(1, :) .* sin(est_traj_B(2, :));

error_x        = x_B_true - x_B_est_cart;
error_y        = y_B_true - y_B_est_cart;
combined_error = sqrt(error_x.^2 + error_y.^2);

% Calculate mean errors
mean_err_x        = mean(error_x);
mean_err_y        = mean(error_y);
mean_err_combined = mean(combined_error);

%% 7. Visualizations and Plots

% Plot 1: 2D Cartesian Trajectory Comparison
figure('Name', '2D Cartesian Trajectory Comparison', 'Position', [180, 180, 850, 650]);
plot(x_A_true, y_A_true, 'r--', 'LineWidth', 2); hold on;
plot(x_B_true, y_B_true, 'g-', 'LineWidth', 2);
plot(x_B_est_cart, y_B_est_cart, 'b:', 'LineWidth', 2);
grid on;
title('2D Cartesian Trajectory Comparison (Polar Dynamic Follower)');
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Drone A (Leader) True', 'Drone B (Follower) True (Dynamic)', 'Drone B (Follower) EKF Estimated', 'Location', 'Best');
axis equal;

% Plot 2: Drone B EKF Position Estimation Error with Mean Error Lines
figure('Name', 'Drone B EKF Position Estimation Error', 'Position', [100, 100, 850, 700]);

subplot(3, 1, 1);
plot(time_vector, error_x, 'LineWidth', 1.5, 'Color', '#D95319'); hold on;
plot([time_vector(1), time_vector(end)], [mean_err_x, mean_err_x], 'k--', 'LineWidth', 1.5);
grid on;
title('Drone B EKF Position Estimation Error (with Mean Error Lines)');
ylabel('X-Axis Error (m)');
legend('X Error', ['Mean X Error: ', sprintf('%.2f m', mean_err_x)], 'Location', 'Best');

subplot(3, 1, 2);
plot(time_vector, error_y, 'LineWidth', 1.5, 'Color', '#EDB120'); hold on;
plot([time_vector(1), time_vector(end)], [mean_err_y, mean_err_y], 'k--', 'LineWidth', 1.5);
grid on;
ylabel('Y-Axis Error (m)');
legend('Y Error', ['Mean Y Error: ', sprintf('%.2f m', mean_err_y)], 'Location', 'Best');

subplot(3, 1, 3);
plot(time_vector, combined_error, 'LineWidth', 1.5, 'Color', '#0072BD'); hold on;
plot([time_vector(1), time_vector(end)], [mean_err_combined, mean_err_combined], 'k--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Combined Error (m)');
legend('Combined Error', ['Mean Combined Error: ', sprintf('%.2f m', mean_err_combined)], 'Location', 'Best');

