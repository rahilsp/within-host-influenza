%% Call IR with only V data

%% Import Healthy Mice Data
clear all
load('Handel_data_complete.mat')
load('Results/MLE_Healthy_V/Healthy_V_MLE.mat')
[t1,d1] = dat_func(dat,1);
dat = [t1,d1];
t = linspace(0,12);

[beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k,T0,V0] = Theta_func_new(Theta_s);
Y0 = [T0,0,0,0,V0,0,0];

%% Evaluate Model
sol_1 =  ode45(@(t,Y) IR_new(t,Y,beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k), t, Y0);  
Yf = deval(sol_1,t)';

%% Plot Results
figure(1)
semilogy(t,Yf(:,5),'Linewidth',1.5);
hold on
scatter(t1,d1,80,'markerfacecolor','b','markeredgecolor','k');
xlabel('Time (Days)','FontSize',22)
ylabel('Viral Titer ','FontSize',22)
% ylim([1e3 1e9])

