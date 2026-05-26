function L = fcn_lins_concordance(X,Y)
% fcn_lins_concordance      lins concordance
%
%   L = fcn_lins_concordance(X) calculates the [P x P] concordance matrix
%       between all columns of an [N x P] matrix of observations.
%
%   L = fcn_lins_concordance(X,Y) calculates [1 x 1] concordance between
%       two column vectors X and Y that both have dimensions [N x 1].
%
%   Inputs:         X,      matrix or column vector of observations
%                   Y,      column vector of observations (optional)
%
%   Outputs:        L,      the Lin's concordance measure
%
%   References:     Lawrence, I., & Lin, K. (1989). A concordance 
%                   correlation coefficient to evaluate reproducibility. 
%                   Biometrics, 255-268.
%
%   Richard Betzel, Indiana University, 2021
%

if nargin == 1
    
    [P,Q] = size(X);
    M = mean(X);
    C = cov(X);
    V = var(X);
    D = V + V';
    E = (M - M').^2;
    L = 2*C./(D + E);
    
elseif nargin == 2
    
    [PX,QX] = size(X);
    [PY,QY] = size(Y);
    P = PX;
    MX = mean(X);
    MY = mean(Y);
    X0 = bsxfun(@minus,X,MX);
    Y0 = bsxfun(@minus,Y,MY);
    CXY = (X0'*Y0)/(P - 1);
    
    VX = sum(X0.^2)/(P - 1);
    VY = sum(Y0.^2)/(P - 1);
    
    DXY = bsxfun(@plus,VY(ones(QX,1),:),VX');
    EXY = bsxfun(@minus,MY(ones(QX,1),:),MX').^2;
    L = 2*CXY./(DXY + EXY);
    
end