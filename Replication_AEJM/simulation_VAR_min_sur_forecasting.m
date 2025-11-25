%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Table C-1

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(12345) % Set seed

saveFolder = 'Figures';

%% Simulate data for a VAR(1) with three variables:
numIter = 100;

hor = 40;
minSampleLength = 80;

% T=500; % Length of the data
% T=180; % Length of the data
T= minSampleLength+numIter+hor; % Length of the data

n=2; % variables 

%A=[0.6 -0.3;0.3 0.4]; % Autoregressive matrix A
A=[0.95 -0.3;0.3 0.4]; % Autoregressive matrix A
C=[0.4;0.5];
B0=[0.6 0;-0.3 0.5]; % Structural matrix B0

e=standardize(randn(T,n)); % Generate structural errors N(0,1)

x(:,1)=zeros(n,1)+C; % Data
%x(:,1)=zeros(n,1); % Data
u(:,1)=zeros(n,1); % Reduced-form errors

for t=2:T
        
        u(:,t)=B0*e(t,:)';
        x(:,t)=C+A*x(:,t-1)+u(:,t);
        
end

x = x';
u = u';



%% Estimate a VAR(1) 



iid = [1,2];

p=1;

% Draws set up:
draws=1000;

ndDisc = 10000;
mcmc_const = 25;

for ii = 1:numIter
    vardata = x(1:minSampleLength+ii-1,:);
    y0_bar = mean(vardata,1);
    res_min.(strcat("it_",num2str(ii)))    = bvarGLP_y0(vardata, p, 'noc',0,'sur', 0, 'mcmc', 1,'hz',hor, 'MCMCconst', mcmc_const, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
    res_min_sur.(strcat("it_",num2str(ii)))    = bvarGLP_y0(vardata, p, 'noc',0,'sur', 1, 'mcmc', 1,'hz',hor, 'MCMCconst', mcmc_const, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
    
end
for ii = 1:numIter
    disp(ii)
    disp(res_min.(strcat("it_",num2str(ii))).mcmc.ACCrate)
    disp(res_min_sur.(strcat("it_",num2str(ii))).mcmc.ACCrate)
end
%% Compute RMSE

med_fc_min = nan(hor,n,numIter);
med_fc_min_sur = nan(hor,n,numIter);
fc_error_min = nan(hor,n,numIter);
fc_error_min_sur = nan(hor,n,numIter);
for ii = 1:numIter
    med_fc_min(:,:,ii) = median(res_min.(strcat("it_",num2str(ii))).mcmc.Dforecast,3);
    med_fc_min_sur(:,:,ii) = median(res_min_sur.(strcat("it_",num2str(ii))).mcmc.Dforecast,3);

    fc_error_min(:,:,ii) = med_fc_min(:,:,ii)-x(minSampleLength+ii:minSampleLength+ii+hor-1,:);
    fc_error_min_sur(:,:,ii) = med_fc_min_sur(:,:,ii)-x(minSampleLength+ii:minSampleLength+ii+hor-1,:);
end

se_min = fc_error_min.^2;
se_min_sur = fc_error_min_sur.^2;

mse_min = mean(se_min,3);
mse_min_sur = mean(se_min_sur,3);

rmse_min = sqrt(mse_min);
rmse_min_sur = sqrt(mse_min_sur);

% RMSE table
variables = [1,2];
horizons = [1,4,8,20];

header = {'Minnesota: Y1','Minnesota: Y2','Minnesota and SUR: Y1','Minnesota and SUR: Y2'};
table_body = [rmse_min(horizons,1),rmse_min(horizons,2),rmse_min_sur(horizons,1),rmse_min_sur(horizons,2)];
T = array2table(table_body,'VariableNames',header);
T.Properties.RowNames = arrayfun(@num2str, horizons, 'UniformOutput', false);
disp(T)


writetable(T, 'Figures/TabC1_RMSE.csv');

