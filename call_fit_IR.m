load('Handel_data_complete.mat')
load('Results/IR_MLE_All_Data/Pooled_MLE.mat')
close all
%%%%%%%%% Set-Up %%%%%%%%%%%%%%%%%%%
u = 10;
b = 3;

%% Parameter Values

Theta = Theta_s;

%%%%%%%%% Estimate Best-Fit Parameters %%%%%%%%%%%%%%%%%%%

tic
[theta,Theta_s,Y,tf,err,res] = fit_IR_new(dat,Theta,ones(7,1),u,b,1,1);
toc
theta_s = theta;
S1 = ['Error = ',num2str(err)];
disp(S1)
S2 = ['Parameters = '];
disp(S2)
for i = 1:length(theta_s)
    disp(theta_s(i))
end

% save([pwd '/Results/b5/u',num2str(u),'_Results'],'theta','Theta_s','Yf_1','tf_1','Yf_2','tf_2','err','AIC','res')
