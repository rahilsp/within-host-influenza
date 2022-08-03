%% IR Local SA
% z = 1: Pooled
% z = 2: Suppressed
z = 1;
t_end = 30;
N = 100;
sav = 0;

%%
if z == 1
load('Pooled_MLE.mat')
elseif z == 2
load('Suppressed_MLE.mat')
end

load('Handel_data_complete.mat')
[t1,d1,t2,d2,t3,d3,t4,d4,t5,d5,t6,d6,~,~,t8,d8] = dat_func(dat,5);
ti = linspace(0,t_end,500);

%% Set Parameters
[beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k,T0,V0] = Theta_func_new(Theta_s);
tau1 = tau;
s1 = s;
Y0 = [T0,0,0,0,V0,0,0];

s_vec = s1*10.^linspace(-1,1,N);
tau_vec = linspace(0,5,N);

% s_vec = s1;
% tau_vec = 0.5;

l1 = length(s_vec);
l2 = length(tau_vec);
res_mat = zeros(l1,l2,5);

tic
for i = 1:l1
    if floor(i/10) == i/10
        disp(i)
    end
for j = 1:l2
s = s_vec(i);
tau = tau_vec(j);
[t,Y1] = ode45(@(t,Y) IR_new(t,Y,beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k), ti, Y0);
res_mat(j,i,:) = biomarker(t,Y1);
end
end
toc

semilogy(t,Y1(:,5))
max(Y1(:,5))
res_mat(1,1,1);

% vq = griddata(s_vec,tau_vec,res_mat(:,:,1),s1,tau1)

%% SA Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labels = {'Maximum Viral Titer','Time of Maximum Viral Titer','Area Under Viral Titer Curve','Duration of Infection','Fraction of Dead Cells'};
labels1 = {'Max_V','Time_Max_V','AUC','Duration','Frac_Dead'};

for w = 1:5
figure
set(gca,'fontsize',16);
surf(s_vec/s1,tau_vec,res_mat(:,:,w));
xlabel('Relative $s$','Fontsize',26);
ylabel('$\tau$','Fontsize',26);
title(labels(w),'Fontsize',26);
ylim([tau_vec(1), tau_vec(l2)])
colormap jet
shading interp 
view(2)
cb=colorbar;
caxis([min(min((res_mat(:,:,w)))) max(max((res_mat(:,:,w))))])
cb.TickLabelInterpreter='latex';
set(gca, 'YScale', 'linear')
set(gca, 'XScale', 'log')
set(gca,'ColorScale','linear')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
box on
filename = [char(labels1(w)),'_',num2str(z)];
if sav == 1
saveas(gcf,filename,'epsc')
savefig(gcf,filename)
end
% yline(tau1,'k-');
% xline(s1/s1,'k-');
end
shg

% if sav == 1
% filename = ['Heatmap_Innate_N',num2str(N)];
% save(filename)
% end
