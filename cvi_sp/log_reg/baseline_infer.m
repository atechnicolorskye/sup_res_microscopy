function [nlz, log_loss, time] = baseline_infer(method_name, y, X, gamma, y_te, X_te, options)
% Baselines for Logistic regression methods

fprintf('%s\n',method_name);
[N,D] = size(X);

% set default options
[maxItersInfer, lowerBoundTol, display, mc, nSamples, beta decay compute_loss] = myProcessOptions(options, ...
'maxItersInfer', 2000, 'lowerBoundTol',1e-4, 'display', 1, 'mc', 0, 'nSamples', 100, 'beta', 0.24, 'decay', 0, 'compute_loss', 0);

% minfunc options
optMinFunc = struct('display', display,...
'Method', 'lbfgs',...
'DerivativeCheck', 'off',...
'LS', 2,...
'MaxIter', maxItersInfer+1,...
'MaxFunEvals', maxItersInfer+1,...
'TolFun', lowerBoundTol,......
'TolX', lowerBoundTol);
w_map = zeros(D,1);

% same initialization as CVI
init_a1 = 0.01;
init_a2 = 0.005;

switch method_name
case 'PG-exact'
    % initialize
    t_lambda = init_a2*ones(N,1);
    tmp1 = bsxfun(@times,X,t_lambda);%diag(t_lambda)*X
    V = inv(diag(gamma) + X'*tmp1);
    dm1 = init_a1*ones(N,1);
    tm1 =X'*dm1;
    m = V*tm1;
    tm = X*m;
    tv = sum(X.*(V*X')',2);
    sW = 1 ./ sqrt(gamma);

    log_loss = zeros(maxItersInfer+1,1);
    nlz = zeros(maxItersInfer+1,1);
    time = zeros(maxItersInfer+1,1);

    if compute_loss==1
        post_dist.mean = m;
        post_dist.covMat =V;
        iter = 0;
        [pred, log_lik]=get_loss(iter, post_dist, X, y, gamma, X_te, y_te);
        log_loss(iter+1) = pred.log_loss;%log_loss0
        nlz(iter+1) = -log_lik;%nlz0
    end

    t_lambda = 0.05*ones(N,1);%much better
    sW = 1 ./ sqrt(gamma);

    r = 1/(beta+1);
    % iterate
    for iter = 1:maxItersInfer
        tic;
        % gradient of fn
        [fb, gmb, gvb] = E_log_p('bernoulli_logit', y, tm, tv, []);%N
        alpha = -gmb; lambda = -2*gvb;
        % gradient of the lower bound
        b = -gamma.*m - X'*alpha;
        tmp1 = bsxfun(@times,X,r*t_lambda);%diag(r*t_lambda)*X
        L = chol(eye(D)+sW*sW'.* (X'*tmp1));
        m = m + (1-r)*sW.*(L\(L'\ (sW.*b) ));
        v = L'\(bsxfun(@times,X',sW));%diag(sW)*X' == bsxfun(@times,X',sW)
        tv = sum(v.^2, 1)';
        tm = X*m;
        t_lambda = r*t_lambda + (1-r)*lambda;
        time(iter+1) = toc;
        
        if compute_loss==1
            %L = chol(eye(D)+sW*sW'.*(X'*diag(t_lambda)*X));
            tmp1 = bsxfun(@times,X,r*t_lambda);%diag(r*t_lambda)*X
            L = chol(eye(D)+sW*sW'.*(X'*tmp1));
            half_V = L'\(diag(sW)); %V=half_V'*half_V; %%V= ( tm2+diag(gamma) )^{-1}
            V = half_V'*half_V;
            % posterior distribution
            post_dist.mean = m;
            post_dist.covMat = V;
            [pred, log_lik]=get_loss(iter, post_dist, X, y, gamma, X_te, y_te);
            log_loss(iter+1) = pred.log_loss;
            nlz(iter+1) = -log_lik;
        end
    end
    time = cumsum(time);

case 'SnK-alg2'
    % initialize
    nrreg = D;
    guessMean = w_map;
    priorPrec = diag(gamma);
    priorMeanTimesPrec =priorPrec*w_map;
    guessPrec = priorPrec;

    %%note that y \in (0,1)
    gradfun = @(xMean,cholP,S)logistic_link(y,X,xMean,cholP,priorMeanTimesPrec,priorPrec,S);

    iterations=maxItersInfer;
    stepSize = beta;

    [nlz log_loss time] = GaussVarApproxHessian_alg2(compute_loss, guessMean,guessPrec,gradfun,iterations,nSamples,stepSize,X,y,X_te,y_te,init_a1,init_a2);

case 'SnK-FG'
    nrreg = D;
    guessMean = w_map;
    priorPrec = diag(gamma);
    priorMeanTimesPrec =priorPrec*w_map;
    guessPrec = priorPrec;

    iterations=maxItersInfer;
    stepSize = beta;

    if isfield(options,'dataset_name')
        switch options.dataset_name
        case 'covtype_binary_scale'
            % This initial value works better
            [nlz log_loss time] =  GaussVarApproxFactorGradientGLM(compute_loss, priorMeanTimesPrec,priorPrec,X,iterations,nSamples,stepSize,X,y,X_te,y_te,0,1);
        otherwise
            [nlz log_loss time] = GaussVarApproxFactorGradientGLM(compute_loss, priorMeanTimesPrec,priorPrec,X,iterations,nSamples,stepSize,X,y,X_te,y_te,init_a1,init_a2);
        end
    else
        [nlz log_loss time] = GaussVarApproxFactorGradientGLM(compute_loss, priorMeanTimesPrec,priorPrec,X,iterations,nSamples,stepSize,X,y,X_te,y_te,init_a1,init_a2);
    end

case 'PG-small-exact'
    % initialize
    t_lambda = init_a2*ones(N,1);
    V = inv(diag(gamma) + X'*diag(t_lambda)*X);
    dm1 = init_a1*ones(N,1);
    tm1 =X'*dm1;
    m = V*tm1;
    tm = X*m;
    tv = sum(X.*(V*X')',2);
    sW0 = 1 ./ sqrt(gamma);
    W = repmat(1./gamma, 1, N).*(X');

    assert( all(gamma>0) )
    half_XX=repmat(sW0,1,N).*(X');
    XX = half_XX'*half_XX;%N by N % X*diag(gamma)*X'
    vv0 = sum(half_XX.*half_XX, 1);
    vv0 = vv0';%diag( X*diag(1./gamma)*X')


    log_loss = zeros(maxItersInfer+1,1);
    nlz = zeros(maxItersInfer+1,1);
    time = zeros(maxItersInfer+1,1);

    if compute_loss==1
        post_dist.mean = m;
        post_dist.covMat =V;
        iter = 0;
        [pred, log_lik]=get_loss(iter, post_dist, X, y, gamma, X_te, y_te);
        log_loss(iter+1) = pred.log_loss;%log_loss0
        nlz(iter+1) = -log_lik;%nlz0
    end

    t_lambda = 0.05*ones(N,1);%much better

    r = 1/(beta+1);
    tic;
    % iterate
    for iter = 1:maxItersInfer
        % gradient of fn
        [fb, gmb, gvb] = E_log_p('bernoulli_logit', y, tm, tv, []);%N
        alpha = -gmb; lambda = -2*gvb;
        dm2 = r*t_lambda;
        sW = sqrt(abs(dm2)).*sign(dm2);
        L = chol(eye(N)+sW*sW'.*XX);
        %v=L'\(repmat(sW,1,N).*XX);
        v=L'\(bsxfun(@times,XX,sW));
        tv = vv0-sum(v.*v,1)';%diag(XX) - diag((XX * ((diag(1./dm2) + XX) \ XX) ));
        tm_part1=L'\(sW.*tm);
        tm = r*tm + (1-r) .* (  -XX*alpha + v'*(v*alpha) + v'*tm_part1 );
        % using linear regression
        %b = -gamma.*m - X'*alpha;
        %A = diag(gamma) + X'*(diag(r*t_lambda)*X);
        %m = m + (1-r)*(A\b);
        %tv = diag(X*(A\(X')));
        %tm = X*m;
        t_lambda = r*t_lambda + (1-r)*lambda;
        if compute_loss==1
            t1=L'\(repmat(sW,1,D).*W');
            m = r*m + (1-r) .* ( t1'*tm_part1 - W*alpha + t1'*(v*alpha)  );
            L = chol(eye(D)+sW0*sW0'.*(X'*diag(t_lambda)*X));
            half_V = L'\(diag(sW0)); %V=half_V'*half_V; %%V= ( tm2+diag(gamma) )^{-1}
            V = half_V'*half_V;
            % posterior distribution
            post_dist.mean = m;
            post_dist.covMat = V;
            [pred, log_lik]=get_loss(iter, post_dist, X, y, gamma, X_te, y_te);
            log_loss(iter+1) = pred.log_loss;
            nlz(iter+1) = -log_lik;
        end
        time(iter+1) = toc;
    end

case 'Chol'
    t_lambda = init_a2*ones(N,1);
    tmp1 = bsxfun(@times,X,t_lambda);%diag(t_lambda)*X
    V = inv(diag(gamma) + X'*tmp1);
    dm1 = init_a1*ones(N,1);
    tm1 =X'*dm1;
    m = V*tm1;

    global times_cache
    times_cache=zeros(maxItersInfer+1,1);
    log_loss = zeros(maxItersInfer+1,1);
    nlz = zeros(maxItersInfer+1,1);

    if compute_loss==1
        post_dist.mean = m;
        post_dist.covMat =V;
        iter = 0;
        [pred, log_lik]=get_loss(iter, post_dist, X, y, gamma, X_te, y_te);
        nlz(iter+1)=-log_lik;%nlz0
        log_loss(iter+1)=pred.log_loss;%log_loss0
    end

    global m_cache;
    global V_cache;
    global iter_counter;
    global is_cache;
    is_cache = compute_loss;
    iter_counter = 0;
    m_cache = zeros(D,maxItersInfer+1);
    V_cache = zeros(D,D,maxItersInfer+1);

    t_lambda = 0.5*ones(N,1);
    tmp1 = bsxfun(@times,X,t_lambda);%diag(t_lambda)*X
    V = inv(diag(gamma) + X'*tmp1);
    dm1 = init_a1*ones(N,1);
    tm1 =X'*dm1;
    m = V*tm1;

    % initialize
    v0 = [m; packcovariance( triu(chol(V)) )];

    % minfunc
    tic
    [v, f, exitflag, inform] = minFunc(@call_log_reg, v0, optMinFunc, y, X, gamma);

    if iter_counter+1<=length(times_cache)
        times_cache(iter_counter+1) = toc;
    else
        times_cache(end) = toc;
    end
    time = times_cache;

    if compute_loss==1
        for ii=2:(maxItersInfer+1)
            idx = ii;
            if idx>iter_counter
                idx = iter_counter;
            end
            post_dist.mean = m_cache(:,idx);
            post_dist.covMat = V_cache(:,:,idx);
            [pred, log_lik]=get_loss(ii-1, post_dist, X, y, gamma, X_te, y_te);
            nlz(ii)=-log_lik;
            log_loss(ii)=pred.log_loss;
        end
    end



otherwise
    error('do not support')
end

end
