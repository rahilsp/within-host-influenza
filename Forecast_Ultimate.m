%% Make a forecast using MCMC or KNN
tstart = tic; 
clc
clear all
close all
load('Pooled_MLE.mat')
label1 = {'TV','TEIV','IR'};
load('Results/Data_06_Oct.mat') % Load training data 


%% Forecast Parameters %%%%%%%%%%
%%% u = 1: TV  u = 2: TEIV  u = 3: IR
u = 3;
gen = 0;
sav = 1;

d_end_vec = [0,6];                % Forecast made at this time 
d = [0.35,0.35,0.35,0.25];        % Distribution of group
K = 6;                            % Number of training data
t_end = 12;                       % Length of simulation
inc = 1;                          % Time increment of observed data 
sigma1 = zeros(4,1);              % Training data noise
sigma2 = zeros(4,1);              % Testing data noise
sigmaL = 0.3;                     % Likelihood function sigma
dis = 1;                          % Display figures?
alpha = 0.5;                      % Forecast figure transparency  
t = [1:inc:t_end]';

%%%%% MCMC Parameter %%%%%%%%%%
tstop1 = 5*60;                    % Time of each prior simulation
tstop2 = 5*60;                    % Time of forecasting simulation
M1 = 1e3;                         % Max MCMC Training Data Simulations
M2 = 1e3;                         % Max MCMC Test Data Simulations
z1 = 1;                           % Autocorrelation sampling of prior
z2 = 1;                           % Autocorrelation sampling of posterior 
n = n_func(u);                    % Prior Historgram Spacing


%% Generate Training Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gen == 1

if K ~= 0
train_vec = zeros(K,t_end/inc);
for i = 1:K
Theta = pheno_func(Theta_s,3,1,d);
V = gen_data_IR(sigma1,t_end,1,Theta);
train_vec(i,:) = V';
end
end


%% Generate Test Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Theta = pheno_func(Theta_s,3,2,0);
V = gen_data_IR(sigma2,t_end,1,Theta)';
dat = [t,V];

% semilogy(train_vec');
% hold on
% scatter(dat(:,1),dat(:,2),50,'filled','k');

end


%% MCMC Training Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('IR_Unknown_Theta.mat')
sd = 0.1*ones(8,1);
tstart = tic;
prior_vec_1 = MCMC_prior(train_vec,t,sigmaL,sd,tstop1,M1,Theta_s,u,n);
prior_vec = thin_func(prior_vec_1,z1);


%% MCMC Test Data & Plot Forecasts
% for ii = 1:length(d_end_vec)
% d_end = d_end_vec(ii);
% if d_end ~= 0
% [theta_vec,acc_vec,guess,MLE] = MCMC_Ult(dat,M,z,sigmaL,sd,tstop2,u,Theta_s,prior_vec,d_end,n);
% else
% theta_vec = prior_vec;
% end

% M2 = 1e3;
% z2 = 1;
% sav = 1;
tf = linspace(0,t_end);
t = [0:inc:t_end];
for ii = 1:length(d_end_vec)
d_end = d_end_vec(ii); 
if d_end ~= 0
[theta_vec_1,acc_vec,guess,MLE] = MCMC_Ult_New(dat,M2,sigmaL,sd,tstop2,u,Theta_s,prior_vec,d_end,n);
else
theta_vec_1 = prior_vec;
end
theta_vec = thin_func(theta_vec_1,z2);
Y_vec = MCMC_sim(theta_vec,tf,t,Theta_s,u);
figure
plot_forecast_MCMC(theta_vec,train_vec,dat,Y_vec,d_end,dis,Theta_s,u,alpha);
if sav == 1
CL = clock;
filename = ['Forecast_',char(label1(u)),'_d',num2str(d_end),'_',date,'_',num2str(CL(4)),num2str(CL(5))];
saveas(gcf,filename,'epsc')
savefig(gcf,filename)
save(filename)
end
end
toc(tstart)


%% Plot individual simulations 
% figure
% hold on
% plot(Y_vec');
% p3 = scatter(dat(:,1),dat(:,2),50,'filled','k');
% set(gca, 'YScale', 'log')
% axisfunc()
% shg

%% Plot training data 
% figure
% t = [1:inc:t_end]';
% plot(t,train_vec,'linewidth',0.5);
% hold on
% scatter(dat(:,1),dat(:,2),50,'filled','k');
% set(gca, 'YScale', 'log')
% axisfunc()
% shg

%% Load data and plot
% plot_forecast_MCMC(theta_vec,train_vec,dat,Y_vec,d_end,dis,Theta_s,u);


%% n_func
function n = n_func(u)

if u == 1 || u == 2
    n = 32;
elseif u == 3    
    n = 12;
end

end
