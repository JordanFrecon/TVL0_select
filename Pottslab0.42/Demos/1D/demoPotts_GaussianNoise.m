%% Potts-denoising, Gaussian noise


lambda = 1;
sigma  = 0.5;

%%
% Load signal
original = loadPcwConst('sample1');
% Add Gaussian noise with parameter sigma = 0.1
f = original + sigma * randn(size(original));

%%
% Compute $L^2$ Potts estimator
figure(1)
tic
pottsL2 = minL2Potts(f, 2*lambda);
toc
showPotts(f, pottsL2, original, 'L2-Potts')

%%

%%

y = f;
xn = pottsL2;

dec.T       = @(x) transformTau_TV(x,1);%
dec.Tadj    = @(x) transformTau_TV_adj(x,1);%
span        = max(f)-min(f);
tol         = 10^5;

1/2*norm(xn - y)^2 + lambda*length(find(abs(dec.T(xn))>span/tol))



addpath(genpath('/Users/Carmino/Dropbox/GALILEO/Matlab_ExtensionEusipco/minimisationl0'));
niter = 10^4;

tic
[xn] = TVl0_vJF_vNewProxL0(y',lambda,niter,mean(y')*ones(size(y')),tol);
toc


1/2*norm(xn - y')^2 + lambda*length(find(abs(dec.T(xn))>span/tol))


figure(1);clf;
plot(f); hold on;
plot(original,'k','linewidth',2);
plot(pottsL2,'r','linewidth',2)
plot(xn,'m','linewidth',2);
legend('y','x','potts','var');
grid on;