function [gm, gv, fn] = E_log_p_mc(obs_, J, m, v, S, D, img_size, scale, psf, background, mode) 
% For simple Poisson
% E_log_p_mc(y, m, v, S, D)
% Based on Khan's E_log_p_mc code from CVI
% Si Kai Lee

    fn = mvnrnd(zeros(J * 2 * D, 1), eye(J * 2 * D), S)';

    % Sample from q  
    if mode == 'chol'
    % Run Cholesky decomposition to enable multivariate sampling
        L = chol(v,'lower');    
        fn = L * fn;
        d2f = zeros(J * 2 * D, J * 2 * D, S);
    elseif mode == 'diag'
        v = diag(v);
        sqrt_v = sqrt(v);
        fn = bsxfun(@times, sqrt_v, fn);
        d2f = zeros(2, J * D, S);
    end

    fn = bsxfun(@plus, m(:), fn);
    % Cannot change values as produces wrong gradients
    % Keep values within limits
    % fn(fn > 32) = 32;
    % fn(fn < 1) = 1;
    % Ensure that values are in 0.5s
    % fn = ceil(fn) - 0.5;
  
    % Compute MC approximation with STORM
    % Initialize collection
    df = zeros(2, J * D, S);
    xs = 0.5:1:img_size;
    [Ygrid, Xgrid] = meshgrid(xs,xs);
    for i = 1:S
        for j = 1:D
            % Create vector for scale
            scale_vec = scale * ones(J, 1); 
            [fn_ , df_, d2f_] = negLikelihood_3_Hes(fn(2*J*j-2*J+1:2*J*j, i), scale_vec, obs_(:, j), psf, background * ones(1024, 1), 1, Xgrid(:), Ygrid(:));
            df_ = reshape(df_, 2, []);
            df(:, J*j-J+1:J*j, i) = df_;
            if mode == 'chol'
                d2f(2*J*j-2*J+1:2*J*j, 2*J*j-2*J+1:2*J*j, i) = d2f_;
            elseif mode == 'diag'
                d2f(:, j, i) = diag(d2f_);
            end
        end
    end
    
    % Compute expectation
    gm = -mean(df, 3);
    gv = -mean(d2f, 3)/2;
end
