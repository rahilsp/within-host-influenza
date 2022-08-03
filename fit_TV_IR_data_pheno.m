%% Fit TV Model to synthetic data from IR model, with varied phenotypes

close all
load('Pooled_MLE.mat')
sav = 1;

t_end = 12;
A = 1e9; % threshhold
sigma = 0*[0.3,0.05,0.1,0];

K = 1e5;
theta_vec = zeros(K,3,4);
Theta_vec = zeros(K,5,4);
err_vec = zeros(K,4);
Theta0 = [4.8370e-06,1.2294e-04,0.6937,7e9,100];

tic
for v = 1:4
    disp(v);
parfor i = 1:K
Theta3 = pheno_func(Theta_s,v,1);
V_vec = gen_data_IR(sigma,t_end,1,Theta3);
dat = [[1:t_end]',V_vec'];
[theta,Theta1,~,ti,Yf,tf,err] = fit_TV(dat,Theta0,3,0,0);
theta_vec(i,:,v) = theta;
Theta_vec(i,:,v) = Theta1;
err_vec(i,v) = err;
shg
end
end
toc


%% Plot Simulations of Viral Titer

% for v = 1:4
% Theta3 = pheno_func(Theta_s,v,2);
% V_vec = gen_data_new(sigma,t_end,1,Theta3);
% dat = [[1:t_end]',V_vec'];
% subplot(2,2,v);
% [theta,Theta1,~,ti,Yf,tf,err] = fit_TV(dat,Theta0,3,1,1);
% ylim([1e1 1e7])
% yticks([1e1 1e4 1e7])
% title(char(label(v)),'FontSize',18);
% xticks([0, 6, 12])
% xlabel('Time (Days)','FontSize',18)
% ylabel('Viral Titer','FontSize',18)
% shg
% end

%% Barchart SSR

figure
boxplot(err_vec,'Symbol','');
ax = gca;
ax.TickLabelInterpreter = 'latex';
label = {'Baseline','Elderly','Obese','Fit'};
xticklabels(label)
xlabel('Characteristic','Fontsize',22)
ylabel('SSR','Fontsize',22')
ylim([0 30])
shg

%% Scatter Plots


figure
for v = 1:4
scatter(theta_vec(:,1,v),theta_vec(:,2,v),5,'filled')
hold on
end
xlabel('$\beta$','Fontsize',26)
ylabel('$p$','Fontsize',26)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend(label,'FontSize',18);
box on
shg

figure
for v = 1:4
scatter(theta_vec(:,1,v),theta_vec(:,3,v),5,'filled')
hold on
end
xlabel('$\beta$','Fontsize',26)
ylabel('$c$','Fontsize',26)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
legend(label,'FontSize',18);
box on
shg

figure
for v = 1:4
scatter(theta_vec(:,2,v),theta_vec(:,3,v),5,'filled')
hold on
end
xlabel('$p$','Fontsize',26)
ylabel('$c$','Fontsize',26)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
legend(label,'FontSize',18);
box on
shg


%% Histograms

figure
x = [theta_vec(:,1,1),theta_vec(:,1,2),theta_vec(:,1,3),theta_vec(:,1,4)];
[~,edges] = histcounts(log10(x));
for v = 1:4
histogram(theta_vec(:,1,v),10.^edges,'Normalization','probability');
hold on
end
set(gca, 'XScale', 'log')
xlabel('$\beta$','Fontsize',26)
% ylabel('Density','Fontsize',22)
legend(label,'FontSize',18);
box on
shg

figure
x = [theta_vec(:,2,1),theta_vec(:,2,2),theta_vec(:,2,3),theta_vec(:,2,4)];
[~,edges] = histcounts(log10(x));
for v = 1:4
histogram(theta_vec(:,2,v),10.^edges,'Normalization','probability');
hold on
end
set(gca, 'XScale', 'log')
xlabel('$p$','Fontsize',26)
% ylabel('Density','Fontsize',22)
legend(label,'FontSize',18);
box on
shg

figure
x = [theta_vec(:,3,1),theta_vec(:,3,2),theta_vec(:,3,3),theta_vec(:,3,4)];
[~,edges] = histcounts(x);
for v = 1:4
histogram(theta_vec(:,3,v),edges,'Normalization','probability');
hold on
end
xlabel('$c$','Fontsize',26)
% ylabel('Density','Fontsize',22)
legend(label,'FontSize',18);
box on
shg


 %% Histograms Lines
% 
% figure
% x = [theta_vec(:,1,1),theta_vec(:,1,2),theta_vec(:,1,3),theta_vec(:,1,4)];
% [~,edges] = histcounts(log10(x));
% for v = 1:4
% [N,edges] = histcounts(theta_vec(:,1,v), 'Normalization','probability');
% edges = edges(2:end) - (edges(2)-edges(1))/2;
% plot(edges, N/max(N),'linewidth',1.5);
% hold on
% end
% set(gca, 'XScale', 'log')
% xlabel('$\beta$','Fontsize',26)
% ylabel('Density','Fontsize',22)
% legend(label,'FontSize',18);
% box on
% shg
% 
% %%
% figure
% for v = 1:4
% [N,edges] = histcounts(theta_vec(:,2,v), 'Normalization','probability');
% edges = edges(2:end) - (edges(2)-edges(1))/2;
% plot(edges, N,'linewidth',1.5);
% hold on
% end
% set(gca, 'XScale', 'log')
% xlabel('$p$','Fontsize',26)
% ylabel('Density','Fontsize',22)
% legend(label,'FontSize',18);
% box on
% shg
% 
% %%
% figure
% for v = 1:4
% [N,edges] = histcounts(theta_vec(:,3,v), 'Normalization','pdf');
% edges = edges(2:end) - (edges(2)-edges(1))/2;
% plot(edges, N,'linewidth',1.5);
% hold on
% end
% xlabel('$c$','Fontsize',26)
% ylabel('Density','Fontsize',22)
% legend(label,'FontSize',18);
% box on
% shg

