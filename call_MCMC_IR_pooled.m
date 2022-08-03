%% Run IR MCMC Simulations

%% Get Data
load('Handel_data_complete.mat')
load('Pooled_MLE.mat')
[t1,d1,t2,d2,t3,d3,t4,d4,t5,d5,t6,d6,~,~,t8,d8] = dat_func(dat,5);
t_fine = linspace(0,t1(length(t1)));
close all
sav = 1;

%% Define parameters from Theta 
Theta = Theta_s;
[beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k,T0,V0] = Theta_func_new(Theta);
Y0 = [T0,0,0,0,V0,0,0];

%% MCMC Parameters
tstop = 180*60;
M = 1e10;
guess = [beta,p,tau,s,r,f,k,V0];
sigma = 1*[0.3,0.05,0.1,0.2];
sd = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.3];
pri = [0, 1e10; 
       0, 1e10;
       0, 10;
       0, 1e10;
       0, 1e10;
       0, 1e10;
       0, 1e10;
       0, 4e4];
z = 1;
b = 5;
u = 4;

tstart = tic;
[theta_vec,pr_vec,L_vec,acc_vec,MLE] = MCMC_IR_v2(dat,M,Theta,pri,guess,z,sigma,sd,b,u,tstop);
acc_prob = size(theta_vec,1)/length(pr_vec);
toc(tstart)

%% Plot MCMC Simulations
figure
plot(theta_vec(:,1),'r')

%% Plot Histograms
sav = 0
par_name = {'$\beta$','$p$','$\tau$','$s$', '$r$', '$f$', '$k$','$V_0$'};
label1 = {'beta','p','tau', 's', 'r', 'f','k','V0'};

for z = 8
figure
histogram(theta_vec(:,z),'Normalization','probability')
xlabel(par_name(z),'FontSize',32)
ylim([0 0.14])
if sav == 1
saveas(gcf,char(label1(z)),'epsc')
end
end

%%
for z = 8
figure
histogram(aa,'Normalization','probability')
xlabel(par_name(z),'FontSize',32)
ylim([0 0.14])
if sav == 1
saveas(gcf,char(label1(z)),'epsc')
end
end

%%
figure
histogram(theta_vec(:,1),'Normalization','probability')
xlabel('$\beta$','FontSize',32)
saveas(gcf,'beta','epsc')

figure
histogram(theta_vec(:,2),'Normalization','probability')
xlabel('$p$','Fontsize',32)
saveas(gcf,'beta','epsc')

figure
histogram(theta_vec(:,3),'Normalization','probability')
xlabel('$\tau$','Fontsize',32)
saveas(gcf,'beta','epsc')

figure
histogram(theta_vec(:,4),'Normalization','probability')
xlabel('$s$','Fontsize',32)
saveas(gcf,'beta','epsc')

figure
histogram(theta_vec(:,5),'Normalization','probability')
xlabel('$r$','Fontsize',32)
saveas(gcf,'beta','epsc')

figure
histogram(theta_vec(:,6),'Normalization','probability')
xlabel('$f$','Fontsize',32)
saveas(gcf,'beta','epsc')

figure
histogram(theta_vec(:,7),'Normalization','probability')
xlabel('$k$','Fontsize',32)
saveas(gcf,'beta','epsc')

figure
histogram(theta_vec(:,8),'Normalization','probability')
xlabel('$V_0$','Fontsize',32)
saveas(gcf,'beta','epsc')

Par_ci = theta_ci(theta_vec,95);
Med = median(theta_vec);


%% Plot Likelihood and Probability 
M1 = length(pr_vec);
figure
scatter([1:M1],L_vec,20,'filled','r')
title('Likelihood')
figure
scatter([1:M1],pr_vec,20,'filled','k')
title('Probability')


%% Run Simulations
K = size(theta_vec,1);
Vi_vec = zeros(K,100,2);
Dead_vec = zeros(K,100,2);
AB_vec = zeros(K,100,2);

for opt = 1:2
for i = 1:K
beta = theta_vec(i,1);
p = theta_vec(i,2);
tau = theta_vec(i,3);
s = theta_vec(i,4);
r = theta_vec(i,5);
f = theta_vec(i,6);
k = theta_vec(i,7);
V0 = theta_vec(i,8);
Y = Vi_func(beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k,T0,V0,t_fine,opt);
Vi_vec(i,:,opt) = Y(:,5);
Dead_vec(i,:,opt) = Y(:,4)/T0;
AB_vec(i,:,opt) = Y(:,7);
end
end

%% Plot Results %%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
[me,yl,yu] = shade_func(Vi_vec(:,:,2),95);
patch([t_fine fliplr(t_fine)], [yl fliplr(yu)],[1 0 0],'Facecolor',[152, 188, 214]/256)
p1 = plot(t_fine,me,'b','linewidth',1.5);
[me,yl,yu] = shade_func(Vi_vec(:,:,1),95);
patch([t_fine fliplr(t_fine)], [yl fliplr(yu)],[1 0 0],'Facecolor',[214,152,208]/256)
p2 = plot(t_fine,me,'r','linewidth',1.5);
p3 = scatter(t5,d5,50,'filled','b');
p4 = scatter(t1,d1,50,'filled','r');
xlabel('Time (Days)','FontSize',22)
ylabel('Viral Titer','FontSize',22)
legend([p3 p4 p1 p2],{'Suppressed Data','Competent Data','Suppressed Model','Competent Model'},'FontSize',18,'location','northwest');
set(gca, 'YScale', 'log')
% ylim([1e1 1e9])
box on
shg


figure
hold on
[me,yl,yu] = shade_func(Dead_vec(:,:,2),95);
patch([t_fine fliplr(t_fine)], [yl fliplr(yu)],[1 0 0],'Facecolor',[152, 188, 214]/256)
p1 = plot(t_fine,me,'b','linewidth',1.5);
[me,yl,yu] = shade_func(Dead_vec(:,:,1),95);
patch([t_fine fliplr(t_fine)], [yl fliplr(yu)],[1 0 0],'Facecolor',[214,152,208]/256)
p2 = plot(t_fine,me,'r','linewidth',1.5);
p3 = scatter(t6,d6,50,'filled','b');
p4 = scatter(t2,d2,50,'filled','r');
xlabel('Time (Days)','FontSize',22)
ylabel('Fraction of Dead Cells','FontSize',22)
legend([p3 p4 p1 p2],{'Suppressed Data','Competent Data','Suppressed Model','Competent Model'},'FontSize',18,'location','northwest');
box on
shg

%%
[me,yl,yu] = shade_func(AB_vec(:,:,1),95);
yl(1) = eps;
yu(1) = eps;
figure
patch([t_fine fliplr(t_fine)], [yl fliplr(yu)],[1 0 0],'Facecolor',[214,152,208]/256)
hold on
p1 = plot(t_fine,me,'r','linewidth',1.5);
p2 = scatter(t3,d3,50,'filled','r');
xlabel('Time (Days)','FontSize',22)
ylabel('Adaptive IR','FontSize',22)
legend([p2 p1],{'Competent Data','Competent Model'},'FontSize',18,'location','southeast');
set(gca, 'YScale', 'log')
ylim([1 1e2])
xlim([0 12])
yline(1);
box on
shg

%%
figure
hold on
p1 = plot(t_fine,IFN(t_fine),'k','linewidth',1.5);
p2 = scatter(t4,d4,50,'filled','r');
p4 = scatter(t8,d8,50,'filled','b');
xlabel('Time (Days)','FontSize',22)
ylabel('Innate IR','FontSize',22)
legend([p4 p2 p1],{'Suppressed Data','Competent Data','Model',},'FontSize',18);
set(gca, 'YScale', 'log')
ylim([1 1e3])
xlim([0 12])
yline(1);
box on
shg
%%

if sav == 1
filename = ['MCMC_Pooled_t',num2str(tstop)];
save(filename)
end

%% Evaluate Models %%%%%%%%%%%%%%%%%%%%%%
function Y = Vi_func(beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k,T0,V0,t_fine,opt)

if opt == 1
sol = ode45(@(t,Y) IR_new(t,Y,beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k), t_fine, [T0,0,0,0,V0,0,0]);  
elseif opt == 2
sol = ode45(@(t,Y) IR_new(t,Y,beta,p,delta,c,kappa,tau,s,0,0,0,0,0), t_fine, [T0,0,0,0,V0,0,0]);  
end
Y = deval(sol,t_fine)';

end
