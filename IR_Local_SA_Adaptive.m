%% IR Local SA
% z = 1: Pooled
% z = 2: Suppressed
% q = 1: k,r
% q = 2: k,f
% q = 3: r,f
z = 1;
q = 2;
t_end = 35;
N = 50;
sav = 1;


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
Y0 = [T0,0,0,0,V0,0,0];

if q == 1
    theta1 = k;
    theta2 = r;
elseif q == 2
    theta1 = k;
    theta2 = f;
elseif q == 3
    theta1 = r;
    theta2 = f;
end

theta1_vec = theta1*10.^linspace(-1,1,N);
theta2_vec = theta2*10.^linspace(-1,1,N);

l1 = length(theta1_vec);
l2 = length(theta2_vec);
res_mat = zeros(l1,l2,5);

tic
W = 1;
for i = 1:l1
    if floor(i/W) == i/W
        disp(i)
    end
for j = 1:l2
    
if q == 1
k = theta1_vec(i);
r = theta2_vec(j);
elseif q == 2
k = theta1_vec(i);
f = theta2_vec(j);
elseif q == 3
r = theta1_vec(i);
f = theta2_vec(j);
end

[t,Y1] = ode45(@(t,Y) IR_new(t,Y,beta,p,delta,c,kappa,tau,s,f,r,m,gamma,k), ti, Y0);
res_mat(j,i,:) = biomarker(t,Y1);
end
end
toc

vq = griddata(theta1_vec,theta2_vec,res_mat(:,:,1),theta1,theta2);

%% SA Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
labels = {'Maximum Viral Titer','Time of Maximum Viral Titer','Area Under Viral Titer Curve','Duration of Infection','Fraction of Dead Cells'};
labels1 = {'Max_V','Time_Max_V','AUC','Duration','Frac_Dead'};

xlab = {'Relative $k$','Relative $k$','Relative $r$'};
ylab = {'Relative $r$','Relative $f$','Relative $f$'};
    
for w = 1:5
figure
set(gca,'fontsize',16);
surf(theta1_vec/theta1,theta2_vec/theta2,res_mat(:,:,w));
xlabel(xlab(q),'Fontsize',26);
ylabel(ylab(q),'Fontsize',26);
title(labels(w),'Fontsize',26);
colormap jet
shading interp 
view(2)
cb=colorbar;
caxis([min(min((res_mat(:,:,w)))) max(max((res_mat(:,:,w))))])
cb.TickLabelInterpreter='latex';
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca,'ColorScale','linear')
box on
filename = [char(labels1(w)),'_q',num2str(q)];
if sav == 1
saveas(gcf,filename,'epsc')
savefig(gcf,filename)
end
end
shg

if sav == 1
filename = ['Heatmap_Adaptive_q',num2str(q),'_N',num2str(N)];
save(filename)
end

