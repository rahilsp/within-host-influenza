%% Run TV MCMC Simulations
% u = 1: TV, u = 2: TEIV, u = 3: IR
clear all
close all
u = 3;
sav = 1;

%% Initalising
load('Handel_data_complete.mat')
load('Healthy_V_MLE.mat')
[t1,d1] = dat_func(dat,1);
dat = [t1,d1];
t = linspace(0,t1(length(t1)));

if u == 1
par_name = {'$\beta$','$p$','$c$','$V_0$'};
par_label = {'beta','p','c','V0'};
elseif u == 2
par_name = {'$\beta$','$p$','$\delta$','$c$','$V_0$'};
par_label = {'beta','p','delta','c','V0'};
elseif u == 3
par_name = {'$\beta$','$p$','$\tau$','$s$','$r$','$f$','$k$','$V_0$'};
par_label = {'beta','p','tau','s','r','f','k','V0'};
end

if u == 3
load('Healthy_V_MLE.mat')
end

%% MCMC Parameters
tstop = 60*60;
M = 1e5;
sigma = 0.3;
sd = 0.1*ones(8,1);

tstart = tic;
[theta_vec,acc_vec,guess,MLE,xx] = MCMC_Ult(dat,M,1,sigma,sd,tstop,u,Theta_s);
toc(tstart)

%% Plot MCMC
plot(acc_vec(:,1))
figure
plot(acc_vec(:,2))

%%
for i = 1:5
figure
plot(theta_vec(:,i))
end

%% Thin Results
% theta_vec_full = theta_vec;
% [theta_vec] = thin_func(theta_vec_full,50);
theta_vec = theta_vec_full;

%% Histograms
close all
sav = 1;
for i = 1:size(theta_vec,2)
figure  
if u == 1 || u == 2 || u == 3
[~,edges] = histcounts(log10(theta_vec(:,i)));
histogram(theta_vec(:,i),10.^edges,'Normalization','probability');
else
histogram(theta_vec(:,i),'Normalization','probability')
end
set(gca,'FontSize',27)
xlabel(par_name(i),'FontSize',45)
ylabel('Normalised Frequency','FontSize',30)
% xlabel(par_name(i),'FontSize',32)
% set(gca, 'XScale', 'log')
if sav == 1
filename = [char(par_label(i)),'_u',num2str(u)];
saveas(gcf,filename,'epsc')
end
end

Par_ci = theta_ci(theta_vec,95);
Med = median(theta_vec);

%% Scatter

% for i = 1:size(theta_vec,2)
% for j = 1:size(theta_vec,2)
%     if i > j
% figure
% scatter(theta_vec(:,j),theta_vec(:,i),7,'filled','b');
% xlabel(par_name(j),'Fontsize',32)
% ylabel(par_name(i),'Fontsize',32)
% box on
% ylabelvert
%     end
% end
% end

%% Run Simulations and calculate DIC
K = size(theta_vec,1);
Y_vec = zeros(K,length(t));
Y_vec2 = zeros(K,length(t1));
D = zeros(K,1);

for i = 1:K
theta = theta_vec(i,:);
[Y,Ya] = ode_func(t,t1,theta,u);
Y_vec(i,:) = Y;
Y_vec2(i,:) = Ya;
D(i) = err_func(Ya,d1,1);
end

DIC2 = mean(D) + 0.5*var(D)

%% Plot
sav = 0;
figure
hold on
[me,yl,yu] = shade_func(Y_vec,95);
patch([t fliplr(t)], [yl fliplr(yu)],[1 0 0],'Facecolor',[152, 188, 214]/256)
p1 = plot(t,me,'b','linewidth',1.5);
p3 = scatter(t1,d1,50,'filled','k');
xlabel('Time (Days)','FontSize',22)
ylabel('Viral Titer','FontSize',22)
box on
legend([p3,p1],{'Data',' Model'},'FontSize',18);
set(gca, 'YScale', 'log')
xlim([0 12])
ylim([1e2 1e6])
yline(1e2,'k')
if sav == 1
filename = ['Vsim_u',num2str(u)];
saveas(gcf,filename,'epsc')
savefig(gcf,filename)
end
shg

% if sav == 1
% filename = ['MCMC_u',num2str(u),'_M',num2str(M)];
% save(filename)
% end


%% ODE Function
function [Y,Ya] = ode_func(t,t1,theta,u)

T0 = 7e9;
kappa = 4;

if u == 1
sol =  ode45(@(t,Y) TV(t,Y,theta(1),theta(2),theta(3)), t, [T0, theta(4)]);     
Y1 = deval(sol,t)';
Y = Y1(:,2);
Y2 = deval(sol,t1)';  
Ya = Y2(:,2);   
elseif u == 2
sol = ode45(@(t,Y) TEIV(t,Y,theta(1),theta(3),theta(2),theta(4),kappa),t, [T0, 0, 0, theta(5)]);
Y1 = deval(sol,t)';
Y = Y1(:,4); 
Y2 = deval(sol,t1)';  
Ya = Y2(:,4);   
elseif u == 3
delta = 2;
c = 10;
m = 1;
gamma = 4.813e-3;     
sol = ode45(@(t,Y) IR_new(t,Y,theta(1),theta(2),delta,c,kappa,theta(3),theta(4),theta(5),theta(6),m,gamma,theta(7)), t, [T0, 0, 0, 0, theta(8), 0, 0]);
Y1 = deval(sol,t)';
Y = Y1(:,5);   
Y2 = deval(sol,t1)';  
Ya = Y2(:,5);   
end
end

