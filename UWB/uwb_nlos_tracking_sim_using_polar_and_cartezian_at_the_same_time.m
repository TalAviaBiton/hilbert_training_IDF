% =========================================================================
% 3D Dynamic Tag Localization: Dual-EKF Hybrid (Cartesian + Polar)
% Features: Runs Cartesian & Polar EKFs in parallel. Dynamically selects
%           the best output using Exponentially Smoothed NIS tracking.
% =========================================================================
clear; clc; close all;

%% 1. Simulation Setup & Configuration
dt = 0.01;                     
T_total = 15;                  
t = 0:dt:T_total;
N = length(t);

%% 2. Ground Truth 3D Motion Models [X; Y; Z]
posA_true = [2*t; 5*sin(0.5*t); 10 + 0.5*t];               
posB_true = [1 + 1.5*t; -2 + 4*cos(0.3*t); 5 + 2*sin(0.2*t)];    

%% 3. Sensor Noise Setup
noise_type  = 'bursty';        
p_burst     = 0.05;            
sigma_burst = 8.0;             

sigma_gps   = 1.2;             
sigma_los   = 0.10;            
mu_nlos     = 1.5;             
sigma_nlos  = 0.40;            

% Generate Noisy Observations
posA_gps = posA_true + sigma_gps * randn(3, N);
true_dist = sqrt(sum((posA_true - posB_true).^2, 1));
[dist_uwb, nlos_states, burst_log] = generate_uwb_channel_3d(...
    true_dist, sigma_los, mu_nlos, sigma_nlos, noise_type, p_burst, sigma_burst);

%% 4. Dual Filter Setup (Cartesian & Polar)
b = 0.96;                      
R_min = sigma_los^2;           
huber_c = 2.0;                 
F = [eye(3), dt*eye(3); zeros(3), eye(3)]; 

% --- Cartesian EKF Initialization ---
x_cart = [0; 0; 0; 0; 0; 0];
P_cart = diag([10, 10, 10, 2, 2, 2]);
Q_cart = diag([0.05, 0.05, 0.05, 0.2, 0.2, 0.2]);
R_cart = R_min;

% --- Polar EKF Initialization ---
pA_0 = posA_gps(:, 1);
rel_0 = posB_true(:, 1) - pA_0; 
r_0 = norm(rel_0);
theta_0 = atan2(rel_0(2), rel_0(1));
phi_0   = atan2(rel_0(3), sqrt(rel_0(1)^2 + rel_0(2)^2));
x_pol = [r_0; theta_0; phi_0; 0; 0; 0];
P_pol = diag([5, 0.5, 0.5, 1, 0.1, 0.1]);
Q_pol = diag([0.02, 0.001, 0.001, 0.1, 0.01, 0.01]);
R_pol = R_min;

% --- Hybrid Selection Parameters ---
alpha_score = 0.15;            % EWMA smoothing factor for filter selection
score_cart = 1.0;              % Initial fitness score
score_pol  = 1.0;

% Logging Allocation
posB_cart_log = zeros(3, N);
posB_pol_log  = zeros(3, N);
posB_best_log = zeros(3, N);   % Final hybrid output
active_filter = zeros(1, N);   % 1 = Cartesian, 2 = Polar

%% 5. Parallel Simulation Loop
for k = 1:N
    pA = posA_gps(:, k);
    z_uwb = dist_uwb(k);
    
    % ==========================================
    % A. Cartesian EKF Processing
    % ==========================================
    x_pred_c = F * x_cart;
    P_pred_c = F * P_cart * F' + Q_cart;
    
    dx = x_pred_c(1) - pA(1); dy = x_pred_c(2) - pA(2); dz = x_pred_c(3) - pA(3);
    dist_pred_c = sqrt(dx^2 + dy^2 + dz^2);
    
    if dist_pred_c > 1e-4
        H_c = [dx/dist_pred_c, dy/dist_pred_c, dz/dist_pred_c, 0, 0, 0];
        y_c = z_uwb - dist_pred_c;
        
        d_k = (1 - b) / (1 - b^k);
        R_samp_c = (y_c^2) - (H_c * P_pred_c * H_c');
        R_cart = max((1 - d_k) * R_cart + d_k * R_samp_c, R_min);
        
        S_nom_c = H_c * P_pred_c * H_c' + R_cart;
        w_c = iff(abs(y_c)/sqrt(S_nom_c) <= huber_c, 1.0, huber_c / (abs(y_c)/sqrt(S_nom_c)));
        S_rob_c = H_c * P_pred_c * H_c' + (R_cart / w_c);
        
        K_c = (P_pred_c * H_c') / S_rob_c;
        x_cart = x_pred_c + K_c * y_c;
        P_cart = (eye(6) - K_c * H_c) * P_pred_c;
        
        nis_cart = (y_c^2) / S_rob_c;
    else
        x_cart = x_pred_c; P_cart = P_pred_c;
        nis_cart = score_cart; % Hold score
    end
    posB_cart_log(:, k) = x_cart(1:3);
    
    % ==========================================
    % B. Polar EKF Processing
    % ==========================================
    x_pred_p = F * x_pol;
    P_pred_p = F * P_pol * F' + Q_pol;
    
    dist_pred_p = x_pred_p(1);
    H_p = [1, 0, 0, 0, 0, 0];
    y_p = z_uwb - dist_pred_p;
    
    if dist_pred_p > 1e-4
        d_k = (1 - b) / (1 - b^k);
        R_samp_p = (y_p^2) - (H_p * P_pred_p * H_p');
        R_pol = max((1 - d_k) * R_pol + d_k * R_samp_p, R_min);
        
        S_nom_p = H_p * P_pred_p * H_p' + R_pol;
        w_p = iff(abs(y_p)/sqrt(S_nom_p) <= huber_c, 1.0, huber_c / (abs(y_p)/sqrt(S_nom_p)));
        S_rob_p = H_p * P_pred_p * H_p' + (R_pol / w_p);
        
        K_p = (P_pred_p * H_p') / S_rob_p;
        x_pol = x_pred_p + K_p * y_p;
        P_pol = (eye(6) - K_p * H_p) * P_pred_p;
        
        nis_pol = (y_p^2) / S_rob_p;
    else
        x_pol = x_pred_p; P_pol = P_pred_p;
        nis_pol = score_pol; % Hold score
    end
    
    % Convert Polar to Cartesian for direct comparison
    r_e = x_pol(1); th_e = x_pol(2); ph_e = x_pol(3);
    posB_pol_cart = pA + [r_e * cos(ph_e) * cos(th_e); 
                          r_e * cos(ph_e) * sin(th_e); 
                          r_e * sin(ph_e)];
    posB_pol_log(:, k) = posB_pol_cart;
    
    % ==========================================
    % C. Filter Selection (Best Outcome Logic)
    % ==========================================
    % Smooth the Normalized Innovation Squared (NIS) to assess filter health
    score_cart = (1 - alpha_score) * score_cart + alpha_score * nis_cart;
    score_pol  = (1 - alpha_score) * score_pol  + alpha_score * nis_pol;
    
    % Select the filter that is currently "healthier" (lower smoothed NIS)
    if score_cart <= score_pol
        posB_best_log(:, k) = posB_cart_log(:, k);
        active_filter(k) = 1;
    else
        posB_best_log(:, k) = posB_pol_log(:, k);
        active_filter(k) = 2;
    end
end

%% 6. Results & Visualization
rmse_cart = sqrt(mean(sum((posB_true - posB_cart_log).^2, 1)));
rmse_pol  = sqrt(mean(sum((posB_true - posB_pol_log).^2, 1)));
rmse_best = sqrt(mean(sum((posB_true - posB_best_log).^2, 1)));

fprintf('===========================================\n');
fprintf('  Cartesian EKF RMSE : %.3f meters\n', rmse_cart);
fprintf('  Polar EKF RMSE     : %.3f meters\n', rmse_pol);
fprintf('  Hybrid Best RMSE   : %.3f meters\n', rmse_best);
fprintf('===========================================\n');

figure('Name', 'Dual-EKF Hybrid Selection', 'Position', [100, 100, 1200, 750]);

% Plot 1: 3D Spatial Trajectories
subplot(2, 2, [1, 3]);
plot3(posA_true(1,:), posA_true(2,:), posA_true(3,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Tag A Ground Truth'); hold on;
plot3(posB_true(1,:), posB_true(2,:), posB_true(3,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Tag B Ground Truth');
plot3(posB_best_log(1,:), posB_best_log(2,:), posB_best_log(3,:), 'r:', 'LineWidth', 2, 'DisplayName', 'Tag B (Hybrid Best)');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('3D Trajectory Tracking (Hybrid Best RMSE: %.2fm)', rmse_best));
legend('Location', 'best'); grid on; view(3);

% Plot 2: Absolute Error Comparison
subplot(2, 2, 2);
err_cart = sqrt(sum((posB_true - posB_cart_log).^2, 1));
err_pol  = sqrt(sum((posB_true - posB_pol_log).^2, 1));
err_best = sqrt(sum((posB_true - posB_best_log).^2, 1));
plot(t, err_cart, 'b-', 'Color', [0 0 1 0.3], 'DisplayName', 'Cartesian Error'); hold on;
plot(t, err_pol, 'r-', 'Color', [1 0 0 0.3], 'DisplayName', 'Polar Error');
plot(t, err_best, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Hybrid Error');
xlabel('Time (s)'); ylabel('Position Error (m)');
title('Instantaneous Error Comparison');
legend('Location', 'best'); grid on;

% Plot 3: Dynamic Filter Selection History
subplot(2, 2, 4);
plot(t, active_filter, 'k-', 'LineWidth', 1.5);
yticks([1 2]); yticklabels({'Cartesian', 'Polar'});
ylim([0.5 2.5]);
xlabel('Time (s)'); ylabel('Active Filter');
title('Real-Time Filter Selection Map');
grid on;

%% Local Functions
function [r_uwb, S_state, burst_log] = generate_uwb_channel_3d(...
    true_dist, sigma_los, mu_nlos, sigma_nlos, noise_type, p_burst, sigma_burst)
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