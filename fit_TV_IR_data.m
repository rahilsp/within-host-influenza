%% Fit TV Model to synthetic data from IR model, with varied immune response parameters 

close all
clear all
for v = 4:5
disp(['v = ', num2str(v)])

label = {'$k$ Scale Factor','$r$ Scale Factor','$f$ Scale Factor','$s$ Scale Factor','$\tau$ (Days)'};
label1 = {'k','r','f','s','tau'};
label2 = {'$k$ SF = ','$r$ SF = ','$f$ SF = ','$s$ SF = ','$\tau$ = '};
sav = 1;
load('Pooled_MLE.mat')

t_end = 12;
A = 1e9; % threshhold
sigma = 0*[0.3,0.05,0.1,0];

inc1 = 1/100;
inc2 = 1/100;
base = Theta_s(ind_func(v));
par_vec_1 = base*10.^[0:inc1:1]';
par_vec_2 = base*10.^[0:-inc2:-1]';
if v == 5
par_vec_1 = base*linspace(1,5/base,50)';
par_vec_2 = base*linspace(1,0,20)';
end
    
l1 = length(par_vec_1);
l2 = length(par_vec_2);
par_vec = [par_vec_1;par_vec_2];
theta_vec = zeros(l1+l2,3);
Theta_vec = zeros(l1+l2,5);
err_vec = zeros(l1+l2,1);
% Theta2 = [1.78e-07,0.0087,10,7e9,100];
Theta2 = [4.8370e-06,1.2294e-04,0.6937,7e9,100];
Theta0 = Theta2;

tic
for i = 1:l1+l2
    if (i/10) == floor(i/10)
        disp(i);
    end
    if i == l1+1    
    Theta0 = Theta2;
    end   
Theta_s(ind_func(v)) = par_vec(i);  
[V_vec, Frac_vec] = gen_data_IR(sigma,t_end,1,Theta_s);
dat = [[1:t_end]',V_vec'];
[theta,Theta1,~,ti,Yf,tf,err] = fit_TV(dat,Theta0,5,0,0,0);
Theta0 = Theta1;
if i == 1
base_T = theta;
end
theta_vec(i,:) = theta;
Theta_vec(i,:) = Theta1;
err_vec(i) = err;
end
toc


[P,E] = sort_func(par_vec,err_vec,l1);
[~,T] = sort_func(par_vec,theta_vec,l1);
[X,Y] = trunc(P,E,T,A);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Simulations of Viral Titer
% 
% par_vec_b = base*[1/10,1/2,2,10];
% if v == 5
% par_vec_b = [0,1,3,5];
% end
%     
% for i = 1:length(par_vec_b)
% Theta_s(ind_func(v)) = par_vec_b(i);  
% [mV,cI] = min(abs(par_vec_b(i)*ones(l1+l2,1)-par_vec));
% Theta0 = Theta_vec(cI,:);
% subplot(2,2,i);
% [V_vec, Frac_vec] = gen_data_IR(sigma,t_end,1,Theta_s);
% dat = [[1:t_end]',V_vec'];
% [theta,Theta1,~,~,~,~,err] = fit_TV(dat,Theta0,5,1,1,0);
% ylim([1e1 1e7])
% yticks([1e1, 1e4, 1e7])
% xticks([0, 6, 12])
% xlabel('Time (Days)','FontSize',18)
% ylabel('Viral Titer','FontSize',18)
% if v ~= 5
% leg = [char(label2(v)),num2str(par_vec_b(i)/base)];
% elseif v == 5
% leg = [char(label2(v)),num2str(par_vec_b(i))];
% end
% title(leg,'FontSize',18);
% if sav == 1
% filename = ['TV_Sim_',char(label1(v))];
% saveas(gcf,filename,'epsc')
% savefig(gcf,filename)
% end
% 
% end

%% Plot Best-Fit Parameters and SSR
% if v == 5
%     base = 1;
% end
% S = good_strip(P/base,E,15);
% 
% figure
% fill([S(1,1) S(1,2) S(1,2) S(1,1)],[eps eps 1e9 1e9],'green', 'FaceAlpha',0.3, 'Edgecolor','none')
% hold on
% p = plot(X/base,Y./base_T,'linewidth',1.5);
% xlabel(label(v),'FontSize',22)
% ylabel('Relative TV Parameter Values','FontSize',22)
% % legend({'$\beta$','$p$','$c$','$V_0$'},'FontSize',22,'Location','Southeast');
% legend(p,{'$\beta$','$p$','$c$'},'FontSize',22,'Location','Southeast');
% % xlim([P(1)/base, P(length(P))/base])
% ylim([min(min(Y./base_T)) max(max(Y./base_T))])
% if v ~= 5
% xlim([10^(-1.1) 10^(1.1)])
% for i = 1:length(par_vec_b)
% xline(par_vec_b(i)/base,'k--','linewidth',1.5);
% end
% set(gca, 'XScale', 'log')
% 
% elseif v == 5
% xlim([-0.1 5.1])
% for i = 1:length(par_vec_b)
% xline(par_vec_b(i),'k--','linewidth',1.5);
% end    
% end
% set(gca, 'YScale', 'log')
% ylim([1e-2 1e2])
% 
% if sav == 1
% filename = ['TV_Par_',char(label1(v))];
% saveas(gcf,filename,'epsc')
% savefig(gcf,filename)
% end


%%
% figure
% fill([S(1,1) S(1,2) S(1,2) S(1,1)],[eps eps 30 30],'green', 'FaceAlpha',0.3, 'Edgecolor','none')
% hold on
% plot(P/base,E,'linewidth',1.5)
% xlabel(label(v),'FontSize',22)
% ylabel('SSR','FontSize',22)
% yline(A,'k--','linewidth',1.5);
% xlim([P(1)/base P(length(P))/base])
% ylim([0 30])
% if v ~= 5
% xlim([10^(-1.1) 10^(1.1)])
%     
% for i = 1:length(par_vec_b)
% xline(par_vec_b(i)/base,'k--','linewidth',1.5);
% end
% set(gca, 'XScale', 'log')
% 
% elseif v == 5
% xlim([-0.1 5.1])
% for i = 1:length(par_vec_b)
% xline(par_vec_b(i),'k--','linewidth',1.5);
% end    
% 
% end

mean_SSR = mean(err_vec);
if sav == 1
filename = ['TV_Err_',char(label1(v))];
% saveas(gcf,filename,'epsc')
% savefig(gcf,filename)
save(filename)
end


end


%%
function ind = ind_func(v)
if v == 1
    ind = 12;
elseif v == 2
    ind = 9;
elseif v == 3
    ind = 8;
elseif v == 4
    ind = 7;
elseif v == 5
    ind = 6;     
end
end