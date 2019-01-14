% demo_TVL0_select_Gaussian.m
%
% Bayesian driven automatic selection of lambda in the 1-dimensional l2-TVl0 minimization problem:
%   min_x  1/2||x-z||_2^2 + lambda||Lx||_0
%   where (Lx_i)_n = x_i(n)-x_i(n-1)
%
% J. Frecon. Version: 20-October-2016.


mydir  = which('demo_TVL0_select_Gaussian.m');
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end));       
addpath(genpath(newdir));
grey = .6.*[1 1 1];


% /!\ Run the following command once per MATLAB session
%installPottslab


%% Parameters

% Synthesis
N           = 1000 ;                    % sample size
SNR         = 1 ;                       % signal to noise ratio
xmax        = 1 ;                       % maximum value of true signal 'x'
xmin        = 0 ;                       % minimum value of true signal 'x'
p           = 0.012;                    % change-point probabilité

% Analysis
LambdaList  = 10.^linspace(-1,1,500);   % range of parameters to examine




%% Synthesis
s           = (xmax-xmin)/(3*SNR);
n           = s*randn(1,N);
K           = p*N;
t           = rand(1,2*K)*2/p;
t           = t;
T           = cumsum(t);
T           = round(T(T<=1.5*N));
mu          = (xmax-xmin)*rand(size(T))+xmin;
x(1:T(1))   = mu(1);
for xx=2:length(T)
    x(T(xx-1)+1:T(xx)) = mu(xx);
end
x           = x(1:N);
y           = x + n;



%% Analysis
[xhat,critlambda,lambdahat] = TVL0_select(y,LambdaList);


%% Display
figure(1);clf;
subplot(211);
plot(y,'color',grey); hold on;
plot(x,'k','linewidth',2);
plot(xhat,'--r','linewidth',2);
l=legend('$y$','$x$','$\widehat{x}$');
set(l,'Interpreter','latex','fontsize',20);
grid on;
subplot(212);
semilogx(LambdaList,critlambda,'r');hold on;
semilogx([lambdahat lambdahat],get(gca,'Ylim'),'--r');
xlabel('$\lambda$','Interpreter','latex');
ylabel('Proposed criterion','Interpreter','latex');
grid on;


