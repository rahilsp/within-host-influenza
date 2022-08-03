%% Fit the TV Model to data using MLE

%% Import Healthy Mice Data
load('Handel_data_complete.mat')
load('Results/MLE_Healthy_V/TV_MLE.mat')
[t1,d1] = dat_func(dat,1);
dat = [t1,d1];

%% Generate IR Data
% load('Results/b5/u4_Results.mat')
% r = Theta_s(7);
% f = 100*Theta_s(8);
% k = Theta_s(9);
% w = Theta_s(13);
% t_end = 10;
% V_vec = gen_data(w,r,k,f,0,0,t_end,1);
% t_dat = 1:t_end;
% dat = [t_dat',V_vec'];

%% Parameters
% beta = 1e-6;
% delta = 2;
% c = 10;
% p = 6e-4;
% T0 = 7e9;
% V0 = 4e4;
% Theta = [beta,p,c,T0,V0];
Theta = Theta_s
q = 1;

hold on
[theta,Theta_s,Y,ti,Yf,tf,err,AIC,res,L] = fit_TV(dat,Theta,q,1,1,1);
err
AIC


