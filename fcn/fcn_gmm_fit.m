function [pi_k, mu_k, Sigma_k, logL, labels, resp, AIC, BIC, logLmax] = fcn_gmm_fit(X, K, maxIter, tol)
% Fast EM algorithm for Gaussian Mixture Models (diagonal covariance)
% with AIC/BIC computed identically to MATLAB's fitgmdist
%
% Inputs:
%   X        [n x d] data matrix
%   K        number of components
%   maxIter  maximum EM iterations (default 100)
%   tol      convergence tolerance on log-likelihood (default 1e-6)
%
% Outputs:
%   pi_k     [1 x K] mixture weights
%   mu_k     [K x d] means
%   Sigma_k  [K x d] diagonal covariance entries
%   logL     log-likelihood trace across iterations
%   labels   [n x 1] hard cluster assignments
%   AIC      Akaike Information Criterion
%   BIC      Bayesian Information Criterion
%   logLmax  exact final log-likelihood

    if nargin < 3, maxIter = 100; end
    if nargin < 4, tol = 1e-6; end

    [n, d] = size(X);

    % ---- Initialization ----
    warning('off', 'stats:kmeans:FailedToConverge'); 
    [labels, mu_k] = kmeans(X, K, 'MaxIter', 100, 'Replicates', 1);
    pi_k = histcounts(labels, 0.5:1:(K+0.5)) / n;
    warning('off', 'stats:kmeans:FailedToConverge'); 

    Sigma_k = zeros(K,d);
    for k = 1:K
        Xk = X(labels == k,:) - mu_k(k,:);
        if ~isempty(Xk)
            Sigma_k(k,:) = var(Xk, 0, 1) + 1e-6;
        else
            Sigma_k(k,:) = ones(1,d); % fallback if empty cluster
        end
    end

    logL = nan(maxIter,1);

    % ---- EM loop ----
    for iter = 1:maxIter
        % --- E-step ---
        logProb = zeros(n,K);
        for k = 1:K
            logProb(:,k) = log(pi_k(k) + eps) + loggausspdf_diag(X, mu_k(k,:), Sigma_k(k,:));
        end

        % Normalize responsibilities
        maxLog = max(logProb,[],2);
        resp = exp(logProb - maxLog);
        resp = resp ./ sum(resp,2);

        % Log-likelihood trace
        logL(iter) = sum(maxLog + log(sum(exp(logProb - maxLog),2)));

        % --- M-step ---
        Nk = sum(resp,1);
        pi_k = Nk / n;
        mu_k = (resp' * X) ./ Nk';

        for k = 1:K
            Xc = X - mu_k(k,:);
            Sigma_k(k,:) = (resp(:,k)' * (Xc.^2)) / Nk(k);
            Sigma_k(k,:) = Sigma_k(k,:) + 1e-6;
        end

        % Convergence check
        if iter > 1 && abs(logL(iter) - logL(iter-1)) < tol
            logL = logL(1:iter);
            break;
        end
    end

    % ---- Final labels ----
    [~, labels] = max(resp, [], 2);

    % ---- Exact log-likelihood ----
    logLmax = computeLogL(X, pi_k, mu_k, Sigma_k);

    % ---- Parameter count (matches MATLAB definition) ----
    numParams = (K - 1) + K*d + K*d; % weights + means + diagonal covariances

    % ---- AIC & BIC ----
    AIC = 2*numParams - 2*logLmax;
    BIC = numParams*log(n) - 2*logLmax;
end

% ============================================================
% Helpers
% ============================================================

function logp = loggausspdf_diag(X, mu, sigma2)
% Log density of multivariate normal with diagonal covariance
    d = size(X,2);
    Xc = X - mu;
    logDet = sum(log(sigma2));
    quad = sum((Xc.^2) ./ sigma2, 2);
    logp = -0.5*(d*log(2*pi) + logDet + quad);
end

function ll = computeLogL(X, pi_k, mu_k, Sigma_k)
% Exact log-likelihood for diagonal GMM
    [n, ~] = size(X);
    K = length(pi_k);
    logProb = zeros(n,K);
    for k = 1:K
        logProb(:,k) = log(pi_k(k) + eps) + loggausspdf_diag(X, mu_k(k,:), Sigma_k(k,:));
    end
    maxLog = max(logProb,[],2);
    ll = sum(maxLog + log(sum(exp(logProb - maxLog),2)));
end
