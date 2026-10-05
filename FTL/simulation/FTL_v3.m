%% =========================================================================
% MULTI-DRONE NAVIGATION AND SIMULATION (LEADER-FOLLOWER)
% =========================================================================
% Description: 
%   - Leader A True Trajectory: Sinusoidal path in 2D Cartesian space.
%   - Configurable Coordinate System: 'polar' or 'cartesian'. When 'polar' is 
%     chosen, all EKF estimation & update calculations run natively in polar.
%   - Visualizations: 
%     1. 2D Cartesian Flight Path View (True Leader, True Follower, Estimate).
%     2. Dedicated Drone B Localization Error Analysis.
% =========================================================================

clear; close all; clc;

%% 1. CONFIGURABLE SIMULATION PARAMETERS
config.coord_system       = 'polar';     % Select: 'polar' or 'cartesian'
config.uwb_representation = 'polar';     % UWB mode: 'polar' or 'cartesian'
config.sep_distance       = 50.0;         % Following distance for Drone B (meters)
config.flight_speed       = 10.0;         % Operational speed (m/s)
config.dt                 = 0.1;         % Time step (seconds)
config.T_total            = 40.0;        % Total simulation time (seconds)

time_vec = 0:config.dt:config.T_total;
num_steps = length(time_vec);

%% 2. TRUE TRAJECTORY GENERATION (SINUSOIDAL LEADER A)
amplitude_y = 6.0;  
freq_y      = 0.15;      

% Define Leader A's True Sinusoidal Path in Cartesian space
true_A_cart.x   = config.flight_speed * time_vec;
true_A_cart.y   = amplitude_y * sin(freq_y * time_vec);
true_A_cart.z   = 5.0 * ones(1, num_steps);
true_A_cart.v   = zeros(1, num_steps);
true_A_cart.psi = zeros(1, num_steps);

for k = 1:num_steps
    if k == 1
        dx = config.flight_speed * config.dt;
        dy = amplitude_y * freq_y * cos(freq_y * time_vec(1)) * config.dt;
    else
        dx = true_A_cart.x(k) - true_A_cart.x(k-1);
        dy = true_A_cart.y(k) - true_A_cart.y(k-1);
    end
    true_A_cart.v(k)   = sqrt(dx^2 + dy^2) / config.dt;
    true_A_cart.psi(k) = atan2(dy, dx);
end

true_A = struct();
true_B = struct();

if strcmpi(config.coord_system, 'polar')
    % --- CONVERT TRUE PATH TO POLAR FOR CALCULATIONS ---
    true_A.r   = sqrt(true_A_cart.x.^2 + true_A_cart.y.^2);
    true_A.th  = atan2(true_A_cart.y, true_A_cart.x);
    true_A.z   = true_A_cart.z;
    true_A.v   = true_A_cart.v;
    true_A.psi = true_A_cart.psi;
    
    for k = 1:num_steps
        lag_steps = round((config.sep_distance / config.flight_speed) / config.dt);
        ref_idx = max(1, k - lag_steps);
        
        xA_ref = true_A_cart.x(ref_idx);
        yA_ref = true_A_cart.y(ref_idx);
        psiA_ref = true_A_cart.psi(ref_idx);
        
        xB_cart = xA_ref - config.sep_distance * cos(psiA_ref);
        yB_cart = yA_ref - config.sep_distance * sin(psiA_ref);
        
        true_B.r(k)   = sqrt(xB_cart^2 + yB_cart^2);
        true_B.th(k)  = atan2(yB_cart, xB_cart);
        true_B.z(k)   = true_A_cart.z(ref_idx);
        true_B.v(k)   = true_A_cart.v(ref_idx);
        true_B.psi(k) = psiA_ref;
    end
else
    % --- CARTESIAN TRUE PATH ---
    true_A = true_A_cart;
    
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
end

%% 3. SENSOR NOISE & MEASUREMENT GENERATION
rng(42); % Seed for reproducibility

gps_noise_std   = 1.0;     
dr_v_noise_std  = 0.2;    
dr_psi_noise_std= 0.05; 
uwb_noise_std   = 0.3;     

meas_A = struct();
meas_B = struct();

if strcmpi(config.coord_system, 'polar')
    meas_A.r     = true_A.r + (gps_noise_std/2) * randn(1, num_steps);
    meas_A.th    = true_A.th + (gps_noise_std / 10.0) * randn(1, num_steps);
    meas_A.baro_z= true_A.z; 
    meas_A.dr_v  = true_A.v + dr_v_noise_std * randn(1, num_steps);
    meas_A.dr_psi= true_A.psi + dr_psi_noise_std * randn(1, num_steps);
    
    true_dist = zeros(1, num_steps);
    for k = 1:num_steps
        xA = true_A.r(k) * cos(true_A.th(k));
        yA = true_A.r(k) * sin(true_A.th(k));
        xB = true_B.r(k) * cos(true_B.th(k));
        yB = true_B.r(k) * sin(true_B.th(k));
        true_dist(k) = sqrt((xA - xB)^2 + (yA - yB)^2 + (true_A.z(k) - true_B.z(k))^2);
    end
    meas_B.uwb_dist     = true_dist + uwb_noise_std * randn(1, num_steps);
    meas_B.baro_z       = true_B.z; 
    meas_B.shared_alt_A = true_A.z; 
    meas_B.dr_v         = true_B.v + dr_v_noise_std * randn(1, num_steps);
    meas_B.dr_psi       = true_B.psi + dr_psi_noise_std * randn(1, num_steps);
else
    meas_A.gps_x  = true_A.x + gps_noise_std * randn(1, num_steps);
    meas_A.gps_y  = true_A.y + gps_noise_std * randn(1, num_steps);
    meas_A.baro_z = true_A.z; 
    meas_A.dr_v   = true_A.v + dr_v_noise_std * randn(1, num_steps);
    meas_A.dr_psi = true_A.psi + dr_psi_noise_std * randn(1, num_steps);
    
    true_dist = sqrt((true_A.x - true_B.x).^2 + (true_A.y - true_B.y).^2 + (true_A.z - true_B.z).^2);
    meas_B.uwb_dist     = true_dist + uwb_noise_std * randn(1, num_steps);
    meas_B.baro_z       = true_B.z; 
    meas_B.shared_alt_A = true_A.z; 
    meas_B.dr_v         = true_B.v + dr_v_noise_std * randn(1, num_steps);
    meas_B.dr_psi       = true_B.psi + dr_psi_noise_std * randn(1, num_steps);
end

%% 4. EXTENDED KALMAN FILTER (EKF) INITIALIZATION
est_A = zeros(5, num_steps);
est_B = zeros(5, num_steps);

if strcmpi(config.coord_system, 'polar')
    est_A(:, 1) = [true_A.r(1); true_A.th(1); true_A.z(1); config.flight_speed; true_A.psi(1)];
    est_B(:, 1) = [true_B.r(1); true_B.th(1); true_B.z(1); config.flight_speed; true_B.psi(1)];
else
    est_A(:, 1) = [true_A.x(1); true_A.y(1); true_A.z(1); config.flight_speed; true_A.psi(1)];
    est_B(:, 1) = [true_B.x(1); true_B.y(1); true_B.z(1); config.flight_speed; true_B.psi(1)];
end

P_A = eye(5) * 1.0;
P_B = eye(5) * 1.0;
P_history_B = zeros(5, num_steps);
P_history_B(:, 1) = diag(P_B);

Q = diag([0.01, 0.01, 0.01, 0.1, 0.05]);
dt = config.dt;

%% 5. EKF ESTIMATION LOOP
for k = 2:num_steps
    if strcmpi(config.coord_system, 'polar')
        %% --- POLAR EKF CALCULATIONS ---
        rA = est_A(1, k-1); thA = est_A(2, k-1); zA = est_A(3, k-1);
        vA = est_A(4, k-1); psiA = est_A(5, k-1);
        
        xA_temp = rA * cos(thA);
        yA_temp = rA * sin(thA);
        xA_pred = xA_temp + vA * cos(psiA) * dt;
        yA_pred = yA_temp + vA * sin(psiA) * dt;
        
        rA_pred  = sqrt(xA_pred^2 + yA_pred^2);
        thA_pred = atan2(yA_pred, xA_pred);
        zA_pred  = zA; 
        vA_pred  = meas_A.dr_v(k);
        psiA_pred= meas_A.dr_psi(k);
        
        x_pred_A = [rA_pred; thA_pred; zA_pred; vA_pred; psiA_pred];
        
        F_A = eye(5); eps_val = 1e-5;
        for j = 1:5
            xp = est_A(:, k-1); xp(j) = xp(j) + eps_val;
            xt = xp(1)*cos(xp(2)) + xp(4)*cos(xp(5))*dt;
            yt = xp(1)*sin(xp(2)) + xp(4)*sin(xp(5))*dt;
            pred_p = [sqrt(xt^2+yt^2); atan2(yt, xt); xp(3); xp(4); xp(5)];
            
            xt0 = rA*cos(thA) + vA*cos(psiA)*dt;
            yt0 = rA*sin(thA) + vA*sin(psiA)*dt;
            pred_0 = [sqrt(xt0^2+yt0^2); atan2(yt0, xt0); zA; vA; psiA];
            F_A(:, j) = (pred_p - pred_0) / eps_val;
        end
        
        P_pred_A = F_A * P_A * F_A' + Q;
        
        z_meas_A = [meas_A.r(k); meas_A.th(k); meas_A.baro_z(k)];
        H_A = [1,0,0,0,0; 0,1,0,0,0; 0,0,1,0,0];
        R_A = diag([(gps_noise_std/2)^2, (gps_noise_std/10)^2, 1e-4]);
        
        y_res_A = z_meas_A - H_A * x_pred_A;
        S_A = H_A * P_pred_A * H_A' + R_A;
        K_A = P_pred_A * H_A' / S_A;
        
        est_A(:, k) = x_pred_A + K_A * y_res_A;
        P_A = (eye(5) - K_A * H_A) * P_pred_A;
        
        % Follower B Polar EKF Update
        rB = est_B(1, k-1); thB = est_B(2, k-1); zB = est_B(3, k-1);
        vB = est_B(4, k-1); psiB = est_B(5, k-1);
        
        xB_temp = rB * cos(thB);
        yB_temp = rB * sin(thB);
        xB_pred = xB_temp + vB * cos(psiB) * dt;
        yB_pred = yB_temp + vB * sin(psiB) * dt;
        
        rB_pred  = sqrt(xB_pred^2 + yB_pred^2);
        thB_pred = atan2(yB_pred, xB_pred);
        zB_pred  = meas_B.shared_alt_A(k); 
        vB_pred  = meas_B.dr_v(k);
        psiB_pred= meas_B.dr_psi(k);
        
        x_pred_B = [rB_pred; thB_pred; zB_pred; vB_pred; psiB_pred];
        
        F_B = eye(5);
        for j = 1:5
            xp = est_B(:, k-1); xp(j) = xp(j) + eps_val;
            xt = xp(1)*cos(xp(2)) + xp(4)*cos(xp(5))*dt;
            yt = xp(1)*sin(xp(2)) + xp(4)*sin(xp(5))*dt;
            pred_p = [sqrt(xt^2+yt^2); atan2(yt, xt); xp(3); xp(4); xp(5)];
            
            xt0 = rB*cos(thB) + vB*cos(psiB)*dt;
            yt0 = rB*sin(thB) + vB*sin(psiB)*dt;
            pred_0 = [sqrt(xt0^2+yt0^2); atan2(yt0, xt0); zB; vB; psiB];
            F_B(:, j) = (pred_p - pred_0) / eps_val;
        end
        
        P_pred_B = F_B * P_B * F_B' + Q;
        
        xA_est = est_A(1, k) * cos(est_A(2, k));
        yA_est = est_A(1, k) * sin(est_A(2, k));
        zA_est = est_A(3, k);
        
        xB_cp = rB_pred * cos(thB_pred);
        yB_cp = rB_pred * sin(thB_pred);
        zB_cp = zB_pred;
        
        dist_est = sqrt((xA_est - xB_cp)^2 + (yA_est - yB_cp)^2 + (zA_est - zB_cp)^2);
        z_meas_B = [meas_B.uwb_dist(k); meas_B.shared_alt_A(k)];
        z_pred_B = [dist_est; zB_cp];
        
        H_B = zeros(2, 5);
        for j = 1:5
            xb_p = x_pred_B; xb_p(j) = xb_p(j) + eps_val;
            xb_c = xb_p(1) * cos(xb_p(2));
            yb_c = xb_p(1) * sin(xb_p(2));
            zb_c = xb_p(3);
            dist_p = sqrt((xA_est - xb_c)^2 + (yA_est - yb_c)^2 + (zA_est - zb_c)^2);
            z_pert = [dist_p; zb_c];
            H_B(:, j) = (z_pert - z_pred_B) / eps_val;
        end
        
        R_B = diag([uwb_noise_std^2, 1e-4]);
        y_res_B = z_meas_B - z_pred_B;
        
        S_B = H_B * P_pred_B * H_B' + R_B;
        K_B = P_pred_B * H_B' / S_B;
        
        est_B(:, k) = x_pred_B + K_B * y_res_B;
        P_B = (eye(5) - K_B * H_B) * P_pred_B;
        P_history_B(:, k) = diag(P_B);
        
    else
        %% --- CARTESIAN EKF CALCULATIONS ---
        xA = est_A(1, k-1); yA = est_A(2, k-1); zA = est_A(3, k-1);
        vA = est_A(4, k-1); psiA = est_A(5, k-1);
        
        xA_pred = xA + vA * cos(psiA) * dt;
        yA_pred = yA + vA * sin(psiA) * dt;
        zA_pred = zA; 
        vA_pred = meas_A.dr_v(k);
        psiA_pred = meas_A.dr_psi(k);
        
        x_pred_A = [xA_pred; yA_pred; zA_pred; vA_pred; psiA_pred];
        
        F_A = [1, 0, 0, cos(psiA)*dt, -vA*sin(psiA)*dt;
               0, 1, 0, sin(psiA)*dt,  vA*cos(psiA)*dt;
               0, 0, 1, 0,            0;
               0, 0, 0, 1,            0;
               0, 0, 0, 0,            1];
           
        P_pred_A = F_A * P_A * F_A' + Q;
        
        z_meas_A = [meas_A.gps_x(k); meas_A.gps_y(k); meas_A.baro_z(k)];
        H_A = [1, 0, 0, 0, 0; 0, 1, 0, 0, 0; 0, 0, 1, 0, 0];
        R_A = diag([gps_noise_std^2, gps_noise_std^2, 1e-4]);
        
        y_res_A = z_meas_A - H_A * x_pred_A;
        S_A = H_A * P_pred_A * H_A' + R_A;
        K_A = P_pred_A * H_A' / S_A;
        
        est_A(:, k) = x_pred_A + K_A * y_res_A;
        P_A = (eye(5) - K_A * H_A) * P_pred_A;

        % Drone B Cartesian EKF Update
        xB = est_B(1, k-1); yB = est_B(2, k-1); zB = est_B(3, k-1);
        vB = est_B(4, k-1); psiB = est_B(5, k-1);
        
        xB_pred = xB + vB * cos(psiB) * dt;
        yB_pred = yB + vB * sin(psiB) * dt;
        zB_pred = meas_B.shared_alt_A(k); 
        vB_pred = meas_B.dr_v(k);
        psiB_pred = meas_B.dr_psi(k);
        
        x_pred_B = [xB_pred; yB_pred; zB_pred; vB_pred; psiB_pred];
        
        F_B = [1, 0, 0, cos(psiB)*dt, -vB*sin(psiB)*dt;
               0, 1, 0, sin(psiB)*dt,  vB*cos(psiB)*dt;
               0, 0, 1, 0,            0;
               0, 0, 0, 1,            0;
               0, 0, 0, 0,            1];
           
        P_pred_B = F_B * P_B * F_B' + Q;
        
        dx_est = est_A(1, k) - xB_pred;
        dy_est = est_A(2, k) - yB_pred;
        dz_est = est_A(3, k) - zB_pred;
        dist_est = sqrt(dx_est^2 + dy_est^2 + dz_est^2);
        
        z_meas_B = [meas_B.uwb_dist(k); meas_B.shared_alt_A(k)];
        z_pred_B = [dist_est; zB_pred];
        
        H_B = [-dx_est/dist_est, -dy_est/dist_est, -dz_est/dist_est, 0, 0;
                0,                0,                1,                  0, 0];
        
        R_B = diag([uwb_noise_std^2, 1e-4]);
        y_res_B = z_meas_B - z_pred_B;
        
        S_B = H_B * P_pred_B * H_B' + R_B;
        K_B = P_pred_B * H_B' / S_B;
        
        est_B(:, k) = x_pred_B + K_B * y_res_B;
        P_B = (eye(5) - K_B * H_B) * P_pred_B;
        P_history_B(:, k) = diag(P_B);
    end
end

%% 6. VISUALIZATION SCRIPTS

% Prepare Cartesian coordinates for true paths and estimation paths
if strcmpi(config.coord_system, 'polar')
    true_A_x = true_A.r .* cos(true_A.th);
    true_A_y = true_A.r .* sin(true_A.th);
    true_B_x = true_B.r .* cos(true_B.th);
    true_B_y = true_B.r .* sin(true_B.th);
    
    est_B_x  = est_B(1,:) .* cos(est_B(2,:));
    est_B_y  = est_B(1,:) .* sin(est_B(2,:));
    est_B_z  = est_B(3,:);
else
    true_A_x = true_A.x;
    true_A_y = true_A.y;
    true_B_x = true_B.x;
    true_B_y = true_B.y;
    
    est_B_x  = est_B(1,:);
    est_B_y  = est_B(2,:);
    est_B_z  = est_B(3,:);
end

% --- Figure 1: 2D Cartesian Flight Path View ---
figure('Name', '2D Cartesian Flight Path View', 'Position', [100, 100, 900, 600], 'Color', 'w');
plot(true_A_x, true_A_y, 'b-', 'LineWidth', 2); hold on;
plot(true_B_x, true_B_y, 'r-', 'LineWidth', 2);
plot(est_B_x,  est_B_y,  'm--', 'LineWidth', 1.5);
grid on; axis equal;
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('True Leader (A)', 'True Follower (B)', 'Follower Estimation (B EKF)', 'Location', 'best');
title(sprintf('2D Cartesian Flight Path (Coord Mode: %s, Sep: %.1fm)', upper(config.coord_system), config.sep_distance));

% --- Figure 2: Dedicated Drone B Localization Error Analysis ---
figure('Name', 'Drone B Localization Error Analysis', 'Position', [150, 150, 1000, 750], 'Color', 'w');

% Calculate physical localization errors in Cartesian metrics
err_xB = true_B_x - est_B_x;
err_yB = true_B_y - est_B_y;
err_zB = true_B.z - est_B_z;
err_pos_2D = sqrt(err_xB.^2 + err_yB.^2);
err_pos_3D = sqrt(err_xB.^2 + err_yB.^2 + err_zB.^2);

subplot(3, 1, 1);
plot(time_vec, err_pos_3D, 'r-', 'LineWidth', 2); hold on;
plot(time_vec, err_pos_2D, 'b--', 'LineWidth', 1.5);
grid on; ylabel('Position Error (m)');
legend('Total 3D Position Error', 'Horizontal (2D) Position Error', 'Location', 'best');
title('Drone B Localization Error (Follower EKF vs True Path)');

subplot(3, 1, 2);
plot(time_vec, err_xB, 'r-', 'LineWidth', 1.2); hold on;
plot(time_vec, err_yB, 'g-', 'LineWidth', 1.2);
grid on; ylabel('Axis Error (m)');
legend('X-Axis Error', 'Y-Axis Error', 'Location', 'best');

subplot(3, 1, 3);
plot(time_vec, err_zB, 'm-', 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('Altitude Error (m)');
legend('Z-Axis (Altitude) Error', 'Location', 'best');