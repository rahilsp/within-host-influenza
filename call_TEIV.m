%% Simulate the TEIV Model

% load('Data/Handel_data_complete.mat')
% load('Results/MLE_Healthy_V/TEIV_MLE.mat')
load('Pooled_MLE.mat')

sigma = zeros(4,1)
t_end = 12;
Theta_s(6) = 0;
V_vec = gen_data_IR(sigma,t_end,1,Theta_s);
% [t1,d1] = dat_func(dat,1);

t1 = 1:12
d1 = V_vec;

%% Parameters
t_end = 13;
t_fine = linspace(0,t_end);

beta = 3e-6;
p = 0.0002;
delta = 0.01;
c = 1;
kappa = 4;
T0 = 7e9;
V0 = 300;
Y0 = [T0,0,0,V0];

[t,Y] = ode45(@(t,Y) TEIV(t,Y,beta,delta,p,c,kappa), t_fine, Y0);

%%%%%%%%%%% Plot Results %%%%%%%%%%%%%%%%%%%%
figure
p1 = semilogy(t,Y(:,4),'linewidth',1);
hold on
set(gca, 'YScale', 'log')
p2 = scatter(t1,d1,50,'k','filled');
xlabel('Time (Days)','FontSize',22)
ylabel('Viral Titer (PFU)','FontSize',22)
legend([p2 p1],{'Data','TEIV Model'},'FontSize',18);
xlim([0 t_end])
shg


