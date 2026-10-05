% =========================================================================
% Monte Carlo Dual-EKF Benchmark (100 Trials: Cartesian vs. Polar EKF)
% Features: 100-Trial Statistical Analysis, Mean/Std/95th Percentile RMSE,
%           Randomized UWB/GPS Noise, Boxplots & CDF Comparison
% =========================================================================
clear; clc; close all;

%% 1. Monte Carlo & Simulation Configuration
mc_trials = 100;               % Number of Monte Carlo iterations
dt = 0.01;                     % Sampling interval (seconds)
T_total = 15;                  % Total simulation duration (seconds)
t = 0:dt:T_total;
N = length(t);

% Ground Truth Trajectories (Deterministic across trials)
posA_true = [2*t; 5*sin(0.5*t); 10 + 0.5*t];               
posB_true = [1 + 1.5*t; -2 + 4*cos(0.3*t); 5 + 2*sin(0.2*t)];    
true_dist = sqrt(sum((posA_true - posB_true).^2, 1));

%% 2. Noise & Kinematic Parameters
noise_type  = 'bursty';        
p_burst     = 0.05;            
sigma_burst = 8.0;             
sigma_gps   = 1.5;             
sigma_los   = 0.10;            
mu_nlos     = 1.5;             
sigma_nlos  = 0.40;            

lim.dist_min_axis = [0.001, 0.001, 0.0005];  
lim.dist_max_axis = [0.080, 0.080, 0.0300];  
lim.vel_min_axis  = [0.050, 0.050, 0.0200];  
lim.vel_max_axis  = [3.000, 3.000, 1.5000];  
lim.acc_min_axis  = [0.000, 0.000, 0.0000];  
lim.acc_max_axis  = [2.000, 2.000, 1.0000];  

lim.dist_min_total = 0.002;                  
lim.dist_max_total = 0.100;                  
lim.vel_min_total  = 0.100;                  
lim.vel_max_total  = 4.000;                  
lim.acc_min_total  = 0.000;                  
lim.acc_max_total  = 2.500;                  

% Preallocate Monte Carlo Metrics Storage
rmse_cart_mc = zeros(mc_trials, 1);
rmse_pol_mc  = zeros(mc_trials, 1);

fprintf('Running Monte Carlo Simulation (%d trials)...\n', mc_trials);

%% 3. Monte Carlo Loop
for mc = 1:mc_trials
    % --- Generate Noise Realizations per Trial ---
    posA_gps_raw = posA_true + sigma_gps * randn(3, N);
    posA_clamped = zeros(3, N);
    velA_clamped = zeros(3, N);
    
    posA_clamped(:, 1) = posA_gps_raw(:, 1);
    velA_clamped(:, 1) = [0; 0; 0];
    for k = 2:N
        [posA_clamped(:, k), velA_clamped(:, k)] = clamp_kinematics(...
            posA_gps_raw(:, k), posA_clamped(:, k-1), velA_clamped(:, k-1), dt, lim);
    end
    
    [dist_uwb, ~, ~] = generate_uwb_channel_3d(...
        true_dist, sigma_los, mu_nlos, sigma_nlos, noise_type, p_burst, sigma_burst);
    
    % --- Filter Initializations ---
    b = 0.96;                      
    R_min = sigma_los^2;           
    huber_c = 2.0;                 
    F = [eye(3), dt*eye(3); zeros(3), eye(3)];
    
    % Cartesian Initialization
    x_cart = [0; 0; 0; 0; 0; 0];
    P_cart = diag([10, 10, 10, 2, 2, 2]);
    Q_cart = diag([0.05, 0.05, 0.05, 0.2, 0.2, 0.2]);
    R_cart = R_min;
    
    % Polar Initialization
    pA_0 = posA_clamped(:, 1);
    rel_0 = posB_true(:, 1) - pA_0; 
    r_0 = norm(rel_0);
    theta_0 = atan2(rel_0(2), rel_0(1));
    phi_0   = atan2(rel_0(3), sqrt(rel_0(1)^2 + rel_0(2)^2));
    
    x_pol = [r_0; theta_0; phi_0; 0; 0; 0];
    P_pol = diag([5, 0.5, 0.5, 1, 0.1, 0.1]);
    Q_pol = diag([0.02, 0.001, 0.001, 0.1, 0.01, 0.01]);
    R_pol = R_min;
    
    posB_cart = zeros(3, N);
    posB_pol  = zeros(3, N);
    
    % --- Single Trial Execution ---
    for k = 1:N
        pA = posA_clamped(:, k);
        z_uwb = dist_uwb(k);
        
        % Cartesian EKF Update
        x_pred_c = F * x_cart;
        P_pred_c = F * P_cart * F' + Q_cart;
        
        dx = x_pred_c(1) - pA(1); dy = x_pred_c(2) - pA(2); dz = x_pred_c(3) - pA(3);
        dist_pred_c = sqrt(dx^2 + dy^2 + dz^2);
        
        if dist_pred_c > 1e-4
            H_c = [dx/dist_pred_c, dy/dist_pred_c, dz/dist_pred_c, 0, 0, 0];
            y_c = z_uwb - dist_pred_c;
            d_k = (1 - b) / (1 - b^k);
            R_sample_c = (y_c^2) - (H_c * P_pred_c * H_c');
            R_cart = max((1 - d_k) * R_cart + d_k * R_sample_c, R_min);
            S_nom_c = H_c * P_pred_c * H_c' + R_cart;
            w_c = iff(abs(y_c)/sqrt(S_nom_c) <= huber_c, 1.0, huber_c / (abs(y_c)/sqrt(S_nom_c)));
            S_rob_c = H_c * P_pred_c * H_c' + (R_cart / w_c);
            K_c = (P_pred_c * H_c') / S_rob_c;
            x_cart = x_pred_c + K_c * y_c;
            P_cart = (eye(6) - K_c * H_c) * P_pred_c;
        else
            x_cart = x_pred_c; P_cart = P_pred_c;
        end
        posB_cart(:, k) = x_cart(1:3);
        
        % Polar EKF Update
        x_pred_p = F * x_pol;
        P_pred_p = F * P_pol * F' + Q_pol;
        
        dist_pred_p = x_pred_p(1);
        H_p = [1, 0, 0, 0, 0, 0];
        y_p = z_uwb - dist_pred_p;
        
        if dist_pred_p > 1e-4
            d_k = (1 - b) / (1 - b^k);
            R_sample_p = (y_p^2) - (H_p * P_pred_p * H_p');
            R_pol = max((1 - d_k) * R_pol + d_k * R_sample_p, R_min);
            S_nom_p = H_p * P_pred_p * H_p' + R_pol;
            w_p = iff(abs(y_p)/sqrt(S_nom_p) <= huber_c, 1.0, huber_c / (abs(y_p)/sqrt(S_nom_p)));
            S_rob_p = H_p * P_pred_p * H_p' + (R_pol / w_p);
            K_p = (P_pred_p * H_p') / S_rob_p;
            x_pol = x_pred_p + K_p * y_p;
            P_pol = (eye(6) - K_p * H_p) * P_pred_p;
        else
            x_pol = x_pred_p; P_pol = P_pred_p;
        end
        
        % Polar to Cartesian Coordinate Mapping
        r_e = x_pol(1); th_e = x_pol(2); ph_e = x_pol(3);
        dx_p = r_e * cos(ph_e) * cos(th_e);
        dy_p = r_e * cos(ph_e) * sin(th_e);
        dz_p = r_e * sin(ph_e);
        posB_pol(:, k) = pA + [dx_p; dy_p; dz_p];
    end
    
    % Store Trial RMSE Values
    err_c = sqrt(sum((posB_true - posB_cart).^2, 1));
    err_p = sqrt(sum((posB_true - posB_pol).^2, 1));
    rmse_cart_mc(mc) = sqrt(mean(err_c.^2));
    rmse_pol_mc(mc)  = sqrt(mean(err_p.^2));
end

%% 4. Statistical Computations & Console Output
mean_cart = mean(rmse_cart_mc);
std_cart  = std(rmse_cart_mc);
p95_cart  = prctile(rmse_cart_mc, 95);

mean_pol = mean(rmse_pol_mc);
std_pol  = std(rmse_pol_mc);
p95_pol  = prctile(rmse_pol_mc, 95);

fprintf('\n=========================================================================\n');
fprintf('       MONTE CARLO STATISTICAL BENCHMARK RESULTS (%d TRIALS)\n', mc_trials);
fprintf('=========================================================================\n');
fprintf(' Metric                 | Cartesian EKF (m)     | Polar EKF (m)\n');
fprintf('------------------------+-----------------------+------------------------\n');
fprintf(' Mean RMSE              | %20.4f  | %20.4f\n', mean_cart, mean_pol);
fprintf(' Std Dev RMSE           | %20.4f  | %20.4f\n', std_cart, std_pol);
fprintf(' 95th Percentile RMSE   | %20.4f  | %20.4f\n', p95_cart, p95_pol);
fprintf('=========================================================================\n');

%% 5. Visualization
figure('Name', sprintf('Monte Carlo Analysis (%d Trials)', mc_trials), 'Position', [100, 100, 1100, 500]);

% Subplot 1: Boxplot Comparison
subplot(1, 2, 1);
boxplot([rmse_cart_mc, rmse_pol_mc], 'Labels', {'Cartesian EKF', 'Polar EKF'});
ylabel('3D RMSE (meters)');
title(sprintf('RMSE Distribution Across %d Trials', mc_trials));
grid on;

% Subplot 2: Cumulative Distribution Function (CDF)
subplot(1, 2, 2);
[f_c, x_c] = ecdf(rmse_cart_mc);
[f_p, x_p] = ecdf(rmse_pol_mc);
plot(x_c, f_c, 'b-', 'LineWidth', 2, 'DisplayName', 'Cartesian EKF'); hold on;
plot(x_p, f_p, 'r-', 'LineWidth', 2, 'DisplayName', 'Polar EKF');
xline(p95_cart, 'b--', 'HandleVisibility', 'off');
xline(p95_pol, 'r--', 'HandleVisibility', 'off');
xlabel('Trial RMSE (meters)'); ylabel('Cumulative Probability');
title('Trial RMSE Cumulative Distribution (CDF)');
legend('Location', 'southeast'); grid on;

%% --- HELPER FUNCTIONS ---
function [pos_out, vel_out] = clamp_kinematics(pos_raw, pos_prev, vel_prev, dt, lim)
    d_raw = pos_raw - pos_prev;
    v_raw = d_raw / dt;
    a_raw = (v_raw - vel_prev) / dt;
    
    a_clamped = clamp_vector(a_raw, lim.acc_min_axis, lim.acc_max_axis, lim.acc_min_total, lim.acc_max_total);
    v_target  = vel_prev + a_clamped * dt;
    v_clamped = clamp_vector(v_target, lim.vel_min_axis, lim.vel_max_axis, lim.vel_min_total, lim.vel_max_total);
    d_target  = v_clamped * dt;
    d_clamped = clamp_vector(d_target, lim.dist_min_axis, lim.dist_max_axis, lim.dist_min_total, lim.dist_max_total);
    
    pos_out = pos_prev + d_clamped;
    vel_out = d_clamped / dt;
end

function vec_out = clamp_vector(vec_in, min_axis, max_axis, min_total, max_total)
    vec_out = vec_in;
    for i = 1:3
        mag = abs(vec_out(i));
        sgn = sign(vec_out(i));
        if sgn == 0, sgn = 1; end
        vec_out(i) = sgn * min(max(mag, min_axis(i)), max_axis(i));
    end
    norm_val = norm(vec_out);
    if norm_val > 1e-9
        vec_out = vec_out * (min(max(norm_val, min_total), max_total) / norm_val);
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