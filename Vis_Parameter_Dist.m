%% Visualise Parameter Noise

close all
beta_log = normrnd(log10(7e-8),0.1,1e4,1);
beta_vec = 10.^(beta_log);
figure
histogram(beta_vec,'BinWidth',1e-9)
xlim([1e-8 20e-8])
title('$\beta$')

%%
p_log = normrnd(log10(0.3),0.1,1e4,1);
p_vec = 10.^(p_log);
figure
histogram(p_vec,'BinWidth',0.01)
xlim([0 1])
title('$p$')

%%
s_log = normrnd(log10(0.05),0.2,1e4,1);
s_vec = 10.^(s_log);
figure
histogram(s_vec,'BinWidth',0.005)
% xlim([0 3])
title('$s$')

%%
tau_log = normrnd(log10(1.2),0.1,1e4,1);
tau_vec = 10.^(tau_log);
figure
histogram(tau_vec,'BinWidth',0.02)
xlim([0 3])
title('$\tau$')

%%
V0_log = normrnd(log10(10),0.3,1e4,1);
V0_vec = 10.^(V0_log);
figure
histogram(V0_vec)
% xlim([1 1e4])
title('$V$')

%%

theta_log = normrnd(log10(10),0.1,1e4,1);
theta_vec = 10.^(theta_log);
figure
histogram(theta_vec,'BinWidth',0.1)
xlim([1 20])
title('$c$')

%%
theta_log = normrnd(log10(0.5),0.1,1e4,1);
theta_vec = 10.^(theta_log);
figure
histogram(theta_vec,'BinWidth',0.05)
title('$\delta$')


%%
theta_log = normrnd(log10(1),0.1,1e4,1);
theta_vec = 10.^(theta_log);
figure
histogram(theta_vec,'BinWidth',0.05)
title('$\nu$')

%%
theta_log = normrnd(log10(10),0.1,1e4,1);
theta_vec = 10.^(theta_log);
figure
histogram(theta_vec,'BinWidth',0.05)
title('$\nu$')

%%
theta_log = normrnd(log10(0.2),0.1,1e4,1);
theta_vec = 10.^(theta_log);
figure
histogram(theta_vec,'BinWidth',0.01)
title('$r$')

%%
theta_log = normrnd(log10(4.3e-6),0.2,1e4,1);
theta_vec = 10.^(theta_log);
figure
histogram(theta_vec)
title('$f$')

%%
theta_log = normrnd(log10(8e-4),0.1,1e4,1);
theta_vec = 10.^(theta_log);
figure
histogram(theta_vec)
title('$w$')

%%
theta_log = normrnd(log10(2.5),0.1,1e4,1);
theta_vec = 10.^(theta_log);
figure
histogram(theta_vec)
title('$\mu$')

%%
theta_log = normrnd(log10(2),0.1,1e4,1);
theta_vec = 10.^(theta_log);
figure
histogram(theta_vec)
title('$\delta$')