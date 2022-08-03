%% Fit the IR Model to data of V using MLE

%% Import Healthy Mice Data
load('Handel_data_complete.mat')
load('Healthy_V_MLE.mat')
[t1,d1] = dat_func(dat,1);
dat = [t1,d1];

%% Generate Synthetic Data
% load('Healthy_V_MLE.mat')
% d = 0.35;
% q = 2;
% t_end = 12;
% sigma = zeros(4,1);
% inc = 1;
% t = [1:inc:t_end]';
% Theta = pheno_func(Theta_s,3,q,d);
% V = gen_data_IR(sigma,t_end,1,Theta)';
% dat = [t,V];

%% Parameter Values
% [beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k,T0,V0] = Theta_func_new(Theta);
% beta = 1e-7;
% p = 0.1;
% tau = 1;
% s = 1;
% f = 1;
% r = 1;
% k = 1;
% V0 = 100;
% 
% Theta_s(1) = beta;
% Theta_s(2) = p;
% Theta_s(6) = tau;
% Theta_s(7) = s;
% Theta_s(8) = f;
% Theta_s(9) = r;
% Theta_s(12) = k;
% Theta_s(14) = V0;

%% Call Fit Function 
[theta,Theta,Y,ti,Yf_1,tf_1,err,AIC,res,L] = fit_IR_V(dat,Theta_s,1,1,1);
err
AIC

