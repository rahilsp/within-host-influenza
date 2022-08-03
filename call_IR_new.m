%% Vary parameters of the best-fit IR Model
% z = 1: Healthy
% z = 2: Suppressed
% z = 3: Pooled
z = 3;
t_end = 12;
close all

sigma = 1*[0.3,0.05,0.1,0.2];
% sigma = ones(4,1);

%%
if z == 1
load('Healthy_MLE.mat')
elseif z == 2
load('Suppressed_MLE.mat')
elseif z == 3
load('Pooled_MLE.mat')
end

load('Handel_data_complete.mat')
[t1,d1,t2,d2,t3,d3,t4,d4,t5,d5,t6,d6,~,~,t8,d8] = dat_func(dat,5);
t_fine = linspace(0,t_end,500);

if z == 2
t1 = t5;
d1 = d5;
t2 = t6;
d2 = d6;
t3 = 0;
end


%% Set Parameters
[beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k,T0,V0] = Theta_func_new(Theta_s);
% s = 0.5*s;
% r = 0.5*r;
% k = 0;

Y0 = [T0,0,0,0,V0,0,0];

%% Evaluate Model
sol_1 =  ode45(@(t,Y) IR_new(t,Y,beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k), t_fine, Y0);  
sol_2 =  ode45(@(t,Y) IR_new(t,Y,beta,p,delta,c,kappa,tau,s,0,0,m,gamma,0), t_fine, Y0);  
Yf = deval(sol_1,t_fine)';
Yf2 = deval(sol_2,t_fine)';
tf = t_fine;

%% Calculate SSR
t1a = [0; t1];
t2a = [0; t2];
t3a = [0; t3];
t4a = [0; t4];
t5a = [0; t5];
t6a = [0; t6];
t8a = [0; t8];

Y1 = deval(sol_1,t1a)';
Y2 = deval(sol_1,t2a)';
Y3 = deval(sol_1,t3a)';

C1 = Y1(:,5);
C2 = Y2(:,4)/T0;
C3 = Y3(:,7);

[err1,res1] = err_func(d1,C1([2:length(C1)]),sigma(1),5);
[err2,res2] = err_func(d2,C2([2:length(C2)]),sigma(2),6);
[err3,res3] = err_func(d3,C3([2:length(C3)]),sigma(3),5);

Y5 = deval(sol_2,t5a)';
Y6 = deval(sol_2,t6a)';

C5 = Y5(:,5);
C6 = Y6(:,4)/T0;

[err5,res5] = err_func(d5,C5([2:length(C5)]),sigma(1),5);
[err6,res6] = err_func(d6,C6([2:length(C6)]),sigma(2),6);

err = err1 + err2 + err3 + err5 + err6


%% Plot Results
figure(1)
semilogy(tf,Yf(:,5),'r','Linewidth',1.5);
hold on
semilogy(tf,Yf2(:,5),'b','Linewidth',1.5);
scatter(t1,d1,80,'markerfacecolor','r','markeredgecolor','k');
scatter(t5,d5,80,'markerfacecolor','b','markeredgecolor','k');
xlabel('Time (Days)','FontSize',22)
ylabel('Viral Titer','FontSize',22)
ylim([1e3 1e9])

figure(2)
plot(tf,Yf(:,4)/T0,'r','Linewidth',1.5);
hold on
plot(tf,Yf2(:,4)/T0,'b','Linewidth',1.5);
scatter(t2,d2,80,'markerfacecolor','r','markeredgecolor','k');
scatter(t6,d6,80,'markerfacecolor','b','markeredgecolor','k');
xlabel('Time (Days)','FontSize',22)
ylabel('Fraction of Dead Cells','FontSize',22)

if z == 1 || z == 3
figure(3)
p2 = semilogy(tf,Yf(:,7),'r','Linewidth',1.5);
hold on
scatter(t3,d3,80,'markerfacecolor','r','markeredgecolor','k');
ylim([1 1e2])
xlim([0 t_end])
xlabel('Time (Days)','FontSize',22)
ylabel('Adaptive IR','FontSize',22)
end

figure(4)
p3 = semilogy(t_fine,IFN(t_fine),'Linewidth',1.5);   
hold on
scatter(t4,d4,80,'markerfacecolor','r','markeredgecolor','k');
scatter(t8,d8,80,'markerfacecolor','b','markeredgecolor','k');
xlabel('Time (Days)','FontSize',22)
ylabel('Innate IR','FontSize',22)
ylim([1 1e3])
xlim([0 t_end])
shg

