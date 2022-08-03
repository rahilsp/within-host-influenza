% clear all
% load('/Users/Rahil 1/Dropbox/DPhil/Matlab Project 2/Results/Pheno_Ke5.mat')
% X = theta_vec;

close all
hold on
scatter3(X(:,1,1),X(:,2,1),X(:,3,1),3,'filled')
scatter3(X(:,1,2),X(:,2,2),X(:,3,2),3,'filled')
scatter3(X(:,1,3),X(:,2,3),X(:,3,3),3,'filled')
scatter3(X(:,1,4),X(:,2,4),X(:,3,4),3,'filled')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
% xlim([0 0.5e-4])
% ylim([0 5e-3])
view(45,45)

xlabel('$\hat{\beta}$','Fontsize',26)
ylabel('$\hat{p}$','Fontsize',26)
zlabel('$\hat{c}$','Fontsize',26)
hZLabel = get(gca,'ZLabel');
set(hZLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
shg
