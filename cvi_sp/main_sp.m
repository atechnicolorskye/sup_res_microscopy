% Based on Khan and Lin 2017
% Si Kai Lee

% Missing:
% functions: get_data, log_loss

%% Setup environment and get data
clear all; 
close all;

seed = 1;
[J, mu_1, C_1, A, Q, beta, n_samples, n_iter, img_size, scale, psf, background, mode] = get_params();
% [y, X, y_t, X_t] = get_data(dataset_name);
% [N, D] = size(X);

% Collect log-loss
log_loss = [];

%%%%%%%%%%%%%%%%%%%%%%%
% Single Particle CVI % 
%%%%%%%%%%%%%%%%%%%%%%%
setSeed(seed);
fprintf('CVI for Single Particle\n');

D = 5;

% Generate data using VB_SPT_simulation
all_data = VB_MPT_Simulation(img_size, scale, psf, D, J, mu_1, C_1, Q, A, background, 'Poisson', []);

% Save images to pass to FALCON
filename = VB_WriteImage(all_data.y, './', 'obs.tif', D);  

% Estimate eta
[G, H, mu, Sigma] = get_prior(J, mu_1, C_1, A, Q, D);

% Obtain FALCON input parameters
Gsigma1 = psf; 
Gsigma2 = 100; 
Gsigma_ratio = 1;  
delta_sigma = 0.15;
EM = 0; 
ADU = 1; 
debug = 0;
speed = 'normal';
dummyFrame = 0; 
dfactor = 1;
offset = 1;

% Run FALCON
[fal_output, avg_img, back_est] = FALCON_CPU_rel2_ChangePatameters_SPT('obs.tif', D, dummyFrame, ADU, offset, EM, Gsigma1, Gsigma2, Gsigma_ratio, speed, debug, dfactor);

% Check number of particles per frame, if less duplicate particle with
% highest ratio as most likely to position of multiple particle
num_particle_frame = histc(fal_output(:, 1), unique(fal_output(:, 1)));
for i = 1:D
    if num_particle_frame(i) < J
        num_add = J - num_particle_frame(i);
        [val, idx] = max(fal_output(J*i-(J-1):J*i-num_add, 5));
        rows_add = repmat(fal_output(J*i-(J-idx), :), num_add, 1);
        fal_output = [fal_output(1:J*i-num_add, :); rows_add; fal_output(J*i-(num_add-1):end, :)];
    end
end
        
% Initialize CVI and collection
% Calculate natural parameters of FALCON estimate
% Have to flip as FALCON outputs y then x
est_pos = fliplr(fal_output(:, 2:3))';
ratio = fal_output(:, 5);
est_Sigma = repmat({eye(2)}, J * D, 1);
est_Sigma = cellfun(@times, est_Sigma, mat2cell(ratio, ones(J * D, 1), 1), 'UniformOutput', false);
est_Sigma = blkdiag(est_Sigma{:});
inv_est_Sigma = inv(est_Sigma);
% Add to prior
G = G + inv_est_Sigma * est_pos(:);
H = H - 0.5 .* inv_est_Sigma;

% Create til_lbda(s)
til_lbda_t_1 = zeros(J * 2 * D, 1);
til_lbda_t_2 = zeros(J * 2 * D);
mu_ = zeros(J * 2 * D, n_iter);

% Reshape simulation data
obs_ = reshape(all_data.y, [], D);

% Create block diagonal mask
% Sigma_m = repmat({ones(2)}, D, 1);
% Sigma_m = blkdiag(Sigma_m{:});
    
% Kalman Smoother from PMTK3
% [msmooth, Vsmooth, loglik, VVsmooth] = kalmanSmoother(y', A, eye(2), Q, zeros(2), mu_1', C_1);


%% Run CVI
for iter = 1:n_iter
    % Calculate mean and variance
    Sigma_t = -0.5 .* inv(H + til_lbda_t_2);
    mu_t = Sigma_t * (G + til_lbda_t_1);
    mu_(:, iter) = mu_t;
    mu_t_ = reshape(mu_t, 2, [])';
    
    if mod(iter, n_iter/10) == 0
        figure; 
        plot(all_data.path(:,2),'r','LineWidth',1); hold on; plot(all_data.path(:,1),'b','LineWidth',1); 
        plot(real(mu_t_(:,1)),'k--','LineWidth',1.5); hold on; plot(real(mu_t_(:,2)),'k:','LineWidth',1.5); 
        legend('x','y','x_{est}','y_{est}'); title('Kalman Filter for Single Particle')    
    end
    
    % Step 3
    [df_dm, df_dv, samples] = E_log_p_mc(obs_, J, mu_t, Sigma_t, n_samples, D, img_size, scale, psf, background, mode);
    if mode == 'chol'
        df_dv_mu_t = zeros(J * 2 * D, 1);
        for i = 1:D
            indices = 2*J*i-2*J+1:2*J*i;
            df_dv_mu_t(indices) = df_dv(indices, indices) * mu_t(indices);
        end
        til_lbda_t_2 = (1 - beta) .* til_lbda_t_2 + beta .* df_dv;
    elseif mode == 'diag'
        df_dv_mu_t = df_dv(:) .* mu_t;
        til_lbda_t_2 = (1 - beta) .* til_lbda_t_2 + beta .* diag(df_dv(:));
    end
    til_lbda_t_1 = (1 - beta) .* til_lbda_t_1 + beta .* (df_dm(:) - (2 * df_dv_mu_t));
    
    % Evaluate log-loss
end

%% Calculate results and plot
% Get CVI estimates
var_s_est = - 0.5 .* inv(H + til_lbda_t_2);
s_est = reshape(var_s_est * (G + til_lbda_t_1), 2, [])';

% Tranpose FALCON estimates
est_pos_ = est_pos';

T = reshape(all_data.path', 2, []);
T = T'

figure; 
plot(T(:,1),'r','LineWidth',1); hold on; plot(T(:,2),'b','LineWidth',1); 
plot(real(s_est(:,1)),'k--','LineWidth',1.5); hold on; plot(real(s_est(:,2)),'k:','LineWidth',1.5); 
legend('x','y','x_{est}','y_{est}'); title('Kalman Filter for Single Particle')

% Plot MSEs
mse_fal = mean((est_pos_ - all_data.path) .^ 2);
mse = mean((s_est - all_data.path) .^ 2);
disp(mse); disp(mse_fal);

%% Auxiliary Functions
% Set parameters to desired values
% To further refine
function [J, mu_1, C_1, A, Q, beta, n_samples, n_iter, img_size, scale, psf, background, mode] = get_params()
    J = 3;
    mu_1 = [16, 16];
    C_1 = [3, 0; 0, 3];
    Q = C_1;
    A = eye(2);
    beta = 0.01; % >=0.01 seem to be the best values to match scale of x and y
    n_samples = 10;
    n_iter = 100;
    img_size = 32;
    scale = 2000; % Good to 700
    psf = 1;
    background = 500;
    mode = 'chol';
end

% Get prior
function [G, H, mu, Sigma] = get_prior(J, mu_1, C_1, A, Q, D)
    inv_Q = inv(Q);
    % K
    Q_ = repmat({inv_Q}, 1, J * (D - 1));
    C_ = repmat({inv(C_1)}, 1, J);
    K = blkdiag(C_{:}, Q_{:});
    % L
    AQA = repmat({A * inv_Q * A}, 1, J * (D-1));
    Z_ = repmat({zeros(2)}, 1, J);
    L = blkdiag(AQA{:}, Z_{:});
    % M
    neg_QA_ = repmat({-inv_Q * A}, 1, J);
    neg_QA = blkdiag(neg_QA_{:});
    M = blktridiag(zeros(2 * J), neg_QA, neg_QA, D);
    H_ = K + L + M; 
    % -1/2 * H
    H = - 0.5 * H_;
    % G
    C_1_mu_1 = repmat((C_1 \ mu_1')', J, 1);
    G_ = vertcat(C_1_mu_1, zeros(J * (D-1), 2));
    G = reshape([G_(:, 1), G_(:, 2)]', J * 2 * D, []);
    % Sigma
    Sigma = inv(H_);
    % mu
    mu = H_ \ G;
end