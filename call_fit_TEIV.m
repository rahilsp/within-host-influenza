%% Fit the TEIV Model to data using MLE

%% Import Healthy Mice Data
load('Handel_data_complete.mat')
[t1,d1] = dat_func(dat,1);
dat = [t1,d1];
load('Results/MLE_Healthy_V/TEIV_MLE.mat')

%%
% close all
% load('Pooled_MLE.mat')
% sigma = zeros(4,1);
% t_end = 12;
% % Theta_s(6) = 0;
% Theta_s(6) = 0.5;
% V_vec = gen_data_IR(sigma,t_end,1,Theta_s)';
% t1 = [1:12]';
% d1 = V_vec;
% dat = [t1,d1];

%% Parameters 
% beta = 2e-6;
% p = 0.005;
% delta = 2;
% c = 10;
% kappa = 4;
% T0 = 7e9;
% V0 = 434;
% Theta = [beta,p,delta,c,kappa,T0,V0];
Theta = Theta_s;
q = 1;

hold on
[theta,Theta_s,Y,ti,Yf_1,tf_1,err,AIC,res,L] = fit_TEIV(dat,Theta,q,1,1,1);
set(gca, 'YScale', 'log')
AIC
err

