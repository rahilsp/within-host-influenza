

% legend({'$\hat\beta$','$p$','$c$','$V_0$'},'FontSize',22,'Location','Southeast');
legend({'$\hat\beta$','$\hat{p}$','$\hat{c}$','$\hat{V_0}$'},'FontSize',22);
% xlabel('$1/\tilde{\delta}$ (Days)')
xlabel('$1/\tilde{\kappa}$ (Days)')



% ylabel('SSR','Fontsize',30)
% xlabel('Characteristic','Fontsize',30)
% 
% %%
% set(gca,'Fontsize',22)
% xlabel('$p$','Fontsize',40)
% ylabel('$c$','Fontsize',40)
% ylabel('Normalised Frequency','FontSize',30)


%%
xticks([1e-7, 1e-6, 1e-5, 1e-4])


% close all
% load('Pooled_MLE.mat')
% t_end = 12;
% sigma = 0*[0.3,0.05,0.1,0];
% sc = 2;
% 
% Theta_s(ind_func(1)) = sc*Theta_s(ind_func(1));
% V_vec = gen_data_IR(sigma,t_end,1,Theta_s);
% plot(V_vec,'b','linewidth',2)
% set(gca, 'YScale', 'log')
% 
% load('Pooled_MLE.mat')
% Theta_s(ind_func(3)) = sc*Theta_s(ind_func(3));
% V_vec = gen_data_IR(sigma,t_end,1,Theta_s);
% hold on
% plot(V_vec,'r--','linewidth',2)
% 
% 
% 
% 
% 
% function ind = ind_func(v)
% if v == 1
%     ind = 12;
% elseif v == 2
%     ind = 9;
% elseif v == 3
%     ind = 8;
% elseif v == 4
%     ind = 7;
% elseif v == 5
%     ind = 6;     
% end
% end