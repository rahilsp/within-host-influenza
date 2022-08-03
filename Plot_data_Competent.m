%% Plot mice data

load('Handel_data_complete.mat')
close all
[t1,d1,t2,d2,t3,d3,t4,d4,t5,d5,t6,d6,~,~,t8,d8] = dat_func(dat,5);

figure
hold on
% scatter(t5,d5,80,'markerfacecolor','b','markeredgecolor','k');
scatter(t1,d1,50,'markerfacecolor','b','markeredgecolor','k');
% plot(t5,d5,'b');
plot(t1,d1,'b');
xlim([0 13])
xlabel('Time (Days)','FontSize',26)
ylabel('Viral Titer (PFU)','FontSize',26)
set(gca, 'YScale', 'log')
ylim([5e2 2e6])
% legend({'Suppressed Data','Competent Data'},'FontSize',18);
box on

% figure
% hold on
% scatter(t6,d6,80,'markerfacecolor','b','markeredgecolor','k');
% scatter(t2,d2,80,'markerfacecolor','r','markeredgecolor','k');
% plot(t6,d6,'b');
% plot(t2,d2,'r');
% xlim([0 12])
% xlabel('Time (Days)','FontSize',22)
% ylabel('Fraction of Dead Cells','FontSize',22)
% legend({'Suppressed Data','Competent Data'},'FontSize',18);
% box on
% 
% figure
% hold on
% scatter(t3,d3,80,'markerfacecolor','r','markeredgecolor','k');
% plot(t3,d3,'r');
% ylim([1 1e2])
% xlim([0 12])
% xlabel('Time (Days)','FontSize',22)
% ylabel('Adaptive IR','FontSize',22)
% set(gca, 'YScale', 'log')
% legend({'Competent Data'},'FontSize',18);
% box on
% 
% figure
% hold on
% scatter(t8,d8,80,'markerfacecolor','b','markeredgecolor','k');
% scatter(t4,d4,80,'markerfacecolor','r','markeredgecolor','k');
% plot(t8,d8,'b');
% plot(t4,d4,'r');
% xlabel('Time (Days)','FontSize',22)
% ylabel('Innate IR','FontSize',22)
% set(gca, 'YScale', 'log')
% ylim([1 1e3])
% xlim([0 12])
% legend({'Suppressed Data','Competent Data'},'FontSize',18);
% box on
