%% Visualise the distriubution of immune parameters for different phenotypes
close all

%% s %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
sc_vec = [1,1,1,1,1,0.5,0.5,2];
alpha = [1,1,1,1,0.9,0.8,0.4,0.9];
mu = 0.0685;
x = linspace(0,0.18,1e3);

col = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980;
       0.9290, 0.6940, 0.1250;
       0.4940, 0.1840, 0.5560;
       0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980;
       0.9290, 0.6940, 0.1250;
       0.4940, 0.1840, 0.5560];
   
   
% col = [0, 0.4470, 0.7410;
%        0.8895, 0.5095, 0.1115;
%        0.8895, 0.5095, 0.1115;
%        0.4940, 0.1840, 0.5560];

for i = 1:8
sc = sc_vec(i);
sigma = 0.12*mu;
y = mygampdf(x,sc*mu,sigma^2);
p(i) = patch([x fliplr(x)], [zeros(1,length(y)) fliplr(y)],col(i,:),'Edgecolor','none','Facealpha',alpha(i));
end

xlabel('$s$','Fontsize',40)
ylabel('Probability Density Function','Fontsize',24)
legend([p(1),p(2),p(3),p(4)],{'Baseline','Elderly','Obese','Partially Immune'},'FontSize',18);
xlim([0 0.18])
ylim([0 75])
box on



%% tau %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
sc_vec = [1,1,1,1,1,1,2,1];
alpha = [1,1,1,1,0.9,0.7,0.9,0.3];
mu = 1.032;
x = linspace(0,3,1e3);

col = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980;
       0.9290, 0.6940, 0.1250;
       0.4940, 0.1840, 0.5560;
       0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980;
       0.9290, 0.6940, 0.1250;
       0.4940, 0.1840, 0.5560];
   
   
% col = [0.4480    0.3187    0.4650;
%        0.4480    0.3187    0.4650;
%        0.9290, 0.6940, 0.1250;
%        0.4940, 0.1840, 0.5560];

for i = 1:8
sc = sc_vec(i);
sigma = 0.12*mu;
y = mygampdf(x,sc*mu,sigma^2);
p(i) = patch([x fliplr(x)], [zeros(1,length(y)) fliplr(y)],col(i,:),'Edgecolor','none','Facealpha',alpha(i));
end

xlabel('$\tau$','Fontsize',40)
ylabel('Probability Density Function','Fontsize',24)
legend([p(1),p(2),p(3),p(4)],{'Baseline','Elderly','Obese','Partially Immune'},'FontSize',18);
ylim([0 5])
xlim([0.5,2.7])
box on
shg

%% k

figure
hold on
sc_vec = [1,1,1,1,1,0.5,1,2];
alpha = [1,1,1,1,0.9,0.9,0.5,0.9];
mu = 2.51;
x = linspace(0,7,1e3);

% col = [0, 0.4470, 0.7410;
%        0.8500, 0.3250, 0.0980;
%        0.9290, 0.6940, 0.1250;
%        0.4940, 0.1840, 0.5560;
%        0, 0.4470, 0.7410;
%        0.8500, 0.3250, 0.0980;
%        0.9290, 0.6940, 0.1250;
%        0.4940, 0.1840, 0.5560];
   
   col = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980;
       0.9290, 0.6940, 0.1250;
       0.4940, 0.1840, 0.5560;
       0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980;
       0.9290, 0.6940, 0.1250;
       0.4940, 0.1840, 0.5560];

% col = [0, 0.4470, 0.7410;
%        0.8500, 0.3250, 0.0980;
%        0.9290, 0.6940, 0.1250;
%        0.4940, 0.1840, 0.5560];

for i = 1:8
sc = sc_vec(i);
sigma = 0.12*mu;
y = mygampdf(x,sc*mu,sigma^2);
p(i) = patch([x fliplr(x)], [zeros(1,length(y)) fliplr(y)],col(i,:),'Edgecolor','none','Facealpha',alpha(i));
end

xlabel('$k$','Fontsize',40)
ylabel('Probability Density Function','Fontsize',24)
legend({'Baseline','Elderly','Obese','Partially Immune'},'FontSize',18);
legend([p(1),p(2),p(3),p(4)],{'Baseline','Elderly','Obese','Partially Immune'},'FontSize',18);
xlim([0 6.5])
ylim([0 2])
box on

% figure
% x = linspace(0,5);
% sc = 1;
% mu = 2.51;
% sigma = 0.1*mu;
% y = mygampdf(x,sc*mu,sigma^2);
% plot(x,y);
% title('$k$')
% 
% hold on
% sc = 0.5;
% mu = 2.51;
% sigma = 0.1*mu;
% y = mygampdf(x,sc*mu,sigma^2);
% plot(x,y);
% title('$k$')
% shg
% 
% hold on
% sc = 2;
% mu = 2.51;
% sigma = 0.1*mu;
% y = mygampdf(x,sc*mu,sigma^2);
% plot(x,y);
% title('$k$')
% shg


