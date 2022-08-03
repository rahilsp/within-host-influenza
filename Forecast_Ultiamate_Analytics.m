%% Plot Paramaters as MCMC chain

X = prior_vec;
% X = theta_vec;

for i = 1:4
    figure
    plot(X(:,i))
    shg
end

%% Histograms

% X = prior_vec;
X = theta_vec;

if u == 1
par_name = {'$\beta$','$p$','$c$','$V_0$'};
par_label = {'beta','p','c','V0'};
elseif u == 2
par_name = {'$\beta$','$p$','$\delta$','$c$','$V_0$'};
par_label = {'beta','p','delta','c','V0'};
elseif u == 3
par_name = {'$\beta$','$p$','$\tau$','$s$','$r$','$f$','$k$','$V_0$'};
par_label = {'beta','p','tau','s','r','f','k','V0'};
end

% X = theta_vec;
sav = 0;
for i = 1:size(X,2)
figure  
if u == 1 || u == 2 || u == 3
[~,edges] = histcounts(log10(X(:,i)));
histogram(X(:,i),10.^edges,'Normalization','probability');
else
histogram(X(:,i),'Normalization','probability')
end
xlabel(par_name(i),'FontSize',32)
set(gca, 'XScale', 'log')
if sav == 1
filename = [char(par_label(i)),'_u',num2str(u)];
saveas(gcf,filename,'epsc')
end
end

Par_ci = theta_ci(X,95);
Med = median(X);

%% Plot Forecast

figure
plot_forecast_MCMC(theta_vec,train_vec,dat,Y_vec,d_end,dis,Theta_s,u,alpha);

%%
figure
plot_forecast_MCMC(prior_vec,train_vec,dat,Y_vec,d_end,dis,Theta_s,u,alpha);

%% Plot Acc Vec
close all
plot(acc_vec(:,1))
figure
plot(acc_vec(:,2))

%% Histogram Prior
n = 128;
[N, edges, mid, loc] = histcn(prior_vec,n,n,n,n);
maxN = max(N,[],'all') 

