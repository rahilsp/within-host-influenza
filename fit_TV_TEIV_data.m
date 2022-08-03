%% Fit TV Model to synthetic data from TEIV model, with varied immune response parameters 

% v = 3: delta
% v = 4: c
% v = 5: kappa
v = 3; 
sav = 1;

close all
load('TEIV_MLE_comp_mice.mat')
label = {'','','$1/\delta$ (Days)','$1/c$ (Days)','$1/\kappa$ (Days)'};
label1 = {'','','delta','c','kappa'};

t_end = 12;
A = 15; % threshhold

inc1 = 0.05;
inc2 = 0.05;
Theta_s = Theta1;
base = Theta_s(v); 
par_vec_1 = base*10.^[0:inc1:0.95]';
par_vec_2 = base*10.^[0:-inc2:-1.5]';

    
l1 = length(par_vec_1);
l2 = length(par_vec_2);
par_vec = [par_vec_1;par_vec_2];
theta_vec = zeros(l1+l2,4);
Theta_vec = zeros(l1+l2,5);
err_vec = zeros(l1+l2,1);
R0_TV = zeros(l1+l2,1);
% Theta2 = [1.78e-07,0.0087,10,7e9,100];
Theta_s = Theta1;
base_T = [Theta1(1),Theta1(2),Theta1(4),Theta1(7)];
Q = Theta1(1)*Theta1(2)*Theta1(6)/Theta1(4);
Theta2 = [3.10e-06,1.7e-4,0.920,7e9,283];
Theta0 = Theta2;


%%
tic
for i = 1:l1+l2
    if (i/10) == floor(i/10)
        disp(i);
    end
    if i == l1+1    
    Theta0 = Theta2;
    end   
Theta_s(v) = par_vec(i);  
V_vec = gen_data_TEIV(0,t_end,1,Theta_s);
dat = [[1:t_end]',V_vec'];
[theta,Theta1,~,~,~,~,err] = fit_TV(dat,Theta0,1,0,0);
Theta0 = Theta1;
theta_vec(i,:) = theta;
Theta_vec(i,:) = Theta1;
err_vec(i) = err;
R0_TV(i) = Theta1(1)*Theta1(2)*Theta1(4)/Theta1(3);
end
toc

[P,E] = sort_func(par_vec,err_vec,l1);
[~,T] = sort_func(par_vec,theta_vec,l1);
[~,R0_TV] = sort_func(par_vec,R0_TV,l1);
[X,Y] = trunc(P,E,T,A);
if v == 3
R0_TEIV = Q./P;
else 
R0_TEIV = Q*ones(length(P),1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
v = 5;
figure
par_vec_2 = [0.25,0.5,1,2];
par_vec_3 = 1./par_vec_2; 
Theta3 = [1.166e-06,0.001016,6.21,7e9,380];
Theta0 = Theta3;
for i = 1:length(par_vec_2)
subplot(2,2,i);
Theta_s(v) = par_vec_3(i);  
V_vec = gen_data_TEIV(0,t_end,1,Theta_s);
dat = [[1:t_end]',V_vec'];
[theta,Theta1,~,~,~,~,err] = fit_TV(dat,Theta0,1,1,1);
if v == 3
title(['$1/\delta$ = ',num2str(par_vec_2(i))],'FontSize',18);
elseif v == 5
title(['$1/\kappa$ = ',num2str(par_vec_2(i))],'FontSize',18);
end   
xticks([0, 6, 12])
xlabel('Time (Days)','FontSize',18)
ylabel('Viral Titer','FontSize',18)
hold on
Theta0 = Theta1;
end
shg

% if v == 3
% legend({['$1/\delta$ = ',num2str(par_vec_2(1))],['$1/\delta$ = ',num2str(par_vec_2(2))],['$1/\delta$ = ',num2str(par_vec_2(3))],['$1/\delta$ = ',num2str(par_vec_2(4))]},'FontSize',18);
% elseif v == 5
% legend({['$1/\kappa$ = ',num2str(par_vec_2(1))],['$1/\kappa$ = ',num2str(par_vec_2(2))],['$1/\kappa$ = ',num2str(par_vec_2(3))],['$1/\delta$ = ',num2str(par_vec_2(4))]},'FontSize',18);
% end

%%

for i = 1:length(par_vec_2)
Theta_s(v) = par_vec_3(i);  
V_vec = gen_data_TEIV(0,t_end,1,Theta_s);
dat = [[1:t_end]',V_vec'];
scatter(dat(:,1),dat(:,2),50,'k','filled');
hold on
end
if sav == 1
filename = ['TV_Sim_TEIV_fit',char(label1(v))];
saveas(gcf,filename,'epsc')
savefig(gcf,filename)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

figure
plot(1./X,Y./base_T,'linewidth',1.5)
hold on
% plot(1./P,R0_TV./R0_TEIV,'linewidth',1.5)
xlabel(label(v),'FontSize',22)
ylabel('Relative TV Parameter Values','FontSize',22)
legend({'$\beta$','$p$','$c$','$V_0$'},'FontSize',22,'Location','Northeast');
% legend({'$\beta$','$p$','$c$'},'FontSize',22,'Location','Southeast');
xlim([0 3])
if v == 3
set(gca, 'YScale', 'log')
end
for i = 1:length(par_vec_2)
xline(par_vec_2(i),'k--','linewidth',1.5);
end
if sav == 1
filename = ['TV_Par_TEIV_fit',char(label1(v))];
saveas(gcf,filename,'epsc')
savefig(gcf,filename)
end

figure
plot(1./P,E,'linewidth',1.5)
xlabel(label(v),'FontSize',22)
ylabel('SSR','FontSize',22)
% yline(A,'k--','linewidth',1.5);
% xlim([P(1)/base P(length(P))/base])
% ylim([0 30])
for i = 1:length(par_vec_2)
xline(par_vec_2(i),'k--','linewidth',1.5);
end
xlim([0 3])
if sav == 1
filename = ['TV_Err_TEIV_fit',char(label1(v))];
saveas(gcf,filename,'epsc')
savefig(gcf,filename)
end

