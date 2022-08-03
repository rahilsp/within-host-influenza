%% Visualise Noise of parameters
load('Data/Handel_data_complete.mat')
sigma = 2*0.3;
[t1,d1,t2,d2,t3,d3,t4,d4,t5,d5,t6,d6,~,~,t8,d8] = dat_func(dat,5);
noise_vec = 1.96*sigma;
yu = 10.^(noise_vec).*d1; 
yl = 10.^(-noise_vec).*d1; 

%%% Plot 
figure
hold on
semilogy(t1,yu,'r')
semilogy(t1,yl,'r')
scatter(t1,d1,50,'k','filled')
set(gca, 'YScale', 'log')
xlabel('Time (Days)','FontSize',22)
ylabel('Viral Titer','FontSize',22)
box on
xlim([0 12])
shg


%% Add Noise
sigma = 0.05;
[t1,d1,t2,d2,t3,d3,t4,d4,t5,d5,t6,d6,~,~,t8,d8] = dat_func(dat,5);
noise_vec = 1.96*sigma;
d = d2;
t = t2;
yu = d+noise_vec; 
yl = d-noise_vec; 

%%% Plot 
figure
hold on
semilogy(t,yu,'r')
semilogy(t,yl,'r')
scatter(t,d,50,'k','filled')
xlabel('Time (Days)','FontSize',22)
ylabel('Viral Titer','FontSize',22)
box on
xlim([0 12])
shg

%% Add Noise
sigma = 0.1;
[t1,d1,t2,d2,t3,d3,t4,d4,t5,d5,t6,d6,~,~,t8,d8] = dat_func(dat,5);
noise_vec = 1.96*sigma;
yu = 10.^(noise_vec).*d3; 
yl = 10.^(-noise_vec).*d3; 

%%% Plot 
figure
hold on
semilogy(t3,yu,'r')
semilogy(t3,yl,'r')
scatter(t3,d3,50,'k','filled')
set(gca, 'YScale', 'log')
xlabel('Time (Days)','FontSize',22)
ylabel('Viral Titer','FontSize',22)
box on
xlim([0 12])
shg



%% Add Noise
sigma = 0.2;
[t1,d1,t2,d2,t3,d3,t4,d4,t5,d5,t6,d6,~,~,t8,d8] = dat_func(dat,5);
lt = length(t4);
noise_vec = 1.96*sigma;
yu = 10.^(noise_vec).*d4; 
yl = 10.^(-noise_vec).*d4; 

%%% Plot 
figure
hold on
semilogy(t4,yu,'r')
semilogy(t4,yl,'r')
scatter(t4,d4,50,'k','filled')
set(gca, 'YScale', 'log')
xlabel('Time (Days)','FontSize',22)
ylabel('Viral Titer','FontSize',22)
box on
xlim([0 12])
shg

