%% Simulate the TV Model


%% Generate IR Data
load('Pooled_MLE.mat')
load('Results/MLE_Healthy_V/TV_MLE.mat')
t_end = 12;
M = 1;

sigma = 0*[0.3,0.05,0.1,0];
[V_vec, Frac_vec] = gen_data_new(sigma,t_end,M,Theta_s);
dat = [[1:t_end]',V_vec'];

Theta = [Theta_s(1),Theta_s(2),Theta_s(4),Theta_s(13),Theta_s(14)];

%% Parameters
beta = 0.09*Theta_s(1);
p = 1.1*Theta_s(2);
c = 0.1*Theta_s(4);
T0 = Theta_s(13);
V0 = Theta_s(14);

%% Plot
t_fine = linspace(0,t_end);
[t,Y] = ode45(@(t,Y) TV(t,Y,beta,p,c), t_fine, [T0,V0]);

% %%%%%%%%%%% Plot Results %%%%%%%%%%%%%%%%%%%%
close all
p1 = semilogy(t,Y(:,2),'b','linewidth',1);
hold on
p2 = scatter(dat(:,1),dat(:,2),50,'k','filled');
xlabel('Time (Days)','FontSize',22)
ylabel('Viral Titer','FontSize',22)
legend([p2 p1],{'Data','TV Model'},'FontSize',18);
xlim([0 t_end])
shg


