function d = fcn_jsd(P, Q)
% sqrt_jsd computes the square root of the Jensen-Shannon Divergence
% between two discrete probability distributions P and Q.
%
% Inputs:
%   - P: First probability distribution (vector, sums to 1)
%   - Q: Second probability distribution (vector, sums to 1)
%
% Output:
%   - d: Square root of the Jensen-Shannon divergence (Jensen-Shannon distance)

    % Normalize P and Q to ensure they are probability distributions
    P = P(:) / sum(P);
    Q = Q(:) / sum(Q);

    % Avoid log(0) issues with small epsilon
    epsilon = 1e-12;
    P = P + epsilon;
    Q = Q + epsilon;

    M = 0.5 * (P + Q);  % The average distribution

    % KL divergence helper function
    KL = @(A,B) sum(A .* log2(A ./ B));

    % Compute Jensen-Shannon divergence
    JSD = 0.5 * KL(P, M) + 0.5 * KL(Q, M);

    % Return the square root (Jensen-Shannon distance)
    d = sqrt(JSD);
end
