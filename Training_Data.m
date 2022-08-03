%% Generate Training Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc
close all
load('Pooled_MLE.mat')

M = 1e4;        % Max MCMC Simulations
d = [0.35,0.35,0.35,0.25];        % Distribution of group
d_end = 0;     % Forecast made at this time
K = 100;         % Number of training data
t_end = 12;     % Length of simulation
inc = 1;        % Time increment of observed data 
sigma1 = zeros(4,1);    % Training data noise
sigma2 = zeros(4,1);     % Testing data noise
sigmaL = 0.3;           % Likelihood function sigma
t = [1:inc:t_end]';

train_vec = zeros(K,t_end/inc);
for i = 1:K
Theta = pheno_func(Theta_s,3,1,d);
V = gen_data_IR(sigma1,t_end,1,Theta);
train_vec(i,:) = V';
end

Theta = pheno_func(Theta_s,3,2,0);
V = gen_data_IR(sigma2,t_end,1,Theta)';
dat = [t,V];


t = [1:inc:t_end]';
plot(t,train_vec,'linewidth',0.5);
hold on
scatter(dat(:,1),dat(:,2),50,'filled','k');
set(gca, 'YScale', 'log')
axisfunc()
shg

