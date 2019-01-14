function [xhat,critlambda,lambdahat]= TVL0_select(z,lambdaList, type, param)

% FUNCTION [xhat,critlambda,lambdahat] = TVL0_select(z,lambdaList, type, param)
%
% Bayesian driven automatic selection of 'lambda' in :
%
%   (1-dimensional l2-TVl0 minimization problem)
%       min_x  1/2||x-z||_2^2 + lambda||Lx||_0             
%       where (Lx)_n = x_{n+1}-x_n
%
%   (1-dimensional l1-TVl0 minimization problem)
%       min_x  ||x-z||_1 + lambda||Lx||_0             
%       where (Lx)_n = x_{n+1}-x_n
%
%
% note 1: the TVl0 problems are solved using the Pottslab toolbox developped
% by M. Storath and available at http://pottslab.de/  
%
% note 2: tested and validated on MATLAB R2015a
%
%
% INPUT
%   - 'z' observation (line vector)
%   - 'lambdaList' list of regularization parameters (positive line vector)
%   - 'type' = 'Gaussian' (default) for the l2-TVl0 problem
%              'Laplacian' for the l1-TVl0 problem
%   - 'param' [optional, see article]
%       param.alpha0 = 1 (default)
%       param.alpha1 = 1 (default)
%       param.pi2s02 = 10^4 (default)
%
% OUTPUT 
%   - 'xhat' solution (line vector)
%   - 'critlambda' bayesian driven criterion to minimize (line vector)
%   - 'lambdahat' value of lambda minimizing 'critlambda' and leading to the solution 'xhat'
%
% REFERENCE
%   J. Frecon, N. Pustelnik, N. Dobigeon, H. Wendt, and P. Abry
%   Bayesian selection for the regularization parameter in TVL0 denoising problems
%   Preprint arXiv:1608.07739
%
% J. Frecon
% Version: 20-October-2016


if nargin == 2
    type = 'Gaussian';
end

if nargin <= 3
    param = struct;
end

if ~isfield(param,'alpha0')
    alpha0 = 1;
end

if ~isfield(param,'alpha1')
    alpha1 = 1;
end

if ~isfield(param,'pi2s02')
    A = 10^4;
end


if size(z,1) ~= 1
    z = z';
end


N       = length(z);
span    = max(z)-min(z);
tol     = 10^5;
dec.T   = @(x) [diff(x),0];

if strcmp(type,'Gaussian')
    
    crit    = @(s,lambda,x,z, K) (1/s^2)*( sum((z-x).^2)/2 + lambda*(K-1)) + N*log(2*pi*s^2)/2 + log(s^2) - lambda/(s^2)*(N+alpha0-2) + (N+alpha0-1)*log(A)/2 + (N+alpha0+alpha1-3)*log(1 + exp(lambda/(s^2) - log(A)/2));
    func    = @(a,b) minL2Potts(a, 2*b);
    sigma   = @(a,b) std(a-b);

else
    
    crit    = @(s,lambda,x,z, K) (1/s)*( sum(abs(z-x)) + lambda*(K-1)) + N*log(2*s) + log(s) - lambda/(s)*(N+alpha0-2) + (N+alpha0-1)*log(A) + (N+alpha0+alpha1-3)*log(1 + exp(lambda/(s) - log(A)));
    func    = @(a,b) minL1Potts(a', b)';
    sigma   = @(a,b) sum(abs(a-b))/N;
    
end
    
    
    
for ll=1:length(lambdaList)
    
    lambda           = lambdaList(ll);
    xlambda{ll}      = func(z, lambda);
    
    ind_rlambda      = find(abs(dec.T(xlambda{ll}))>span/tol);
    Klambda          = length(ind_rlambda)+1;
    slambda          = sigma(z,xlambda{ll});
    critlambda(ll)   = crit(slambda,lambda,xlambda{ll},z,Klambda);
    
end
    

[~,indlambdahat]    = min(critlambda);
lambdahat           = lambdaList(indlambdahat);
xhat                = xlambda{indlambdahat};
