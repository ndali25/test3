%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Replication files for paper 
% "Conditional forecasts and scenario analysis with vector autoregressions
% for large cross-sections"
% by M. Banbura, D. Giannone and M. Lenza
% International Journal of Forecasting

% This script calculates the forecasts for 1997-2012 conditional on paths 
% for real GDP, HICP and short-term interest rate
% Parameters of the models are estimated over 1995-2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

warning('off','MATLAB:nearlySingularMatrix')
addpath('../GLPreplicationWeb')
addpath('../GLPreplicationWeb/subroutines')
addpath('../AdditionalFunctions')
addpath('functionsOther')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conditioning variables 
CV = {'YER','HICP','STN'}; % real GDP, HICP and short-term interest rate

%% Start and end of the sample
StartE = [1995 3];
EndE = [2012 12];

%% Number of out-of-sample and initial observations
nQoos = 20; % number of out-of-sample periods at the end of the sample (not included in the estimation sample)
nQinit = 8; % number of initial observations

%% Number of draws
nDraws = 20000; % draws of the parameters

%% Number of lags in the factor VAR and BVAR in differences (in BVAR in
% levels it is p+1)
p = 4; 

%% Settings for the factor model
r = 3; % number of static factors
q = 3; % number of dynamic factors

%% Settings for the BVARs
% Giannone, Lenza, Primiceri procedure settings
Vc = 10e6;    % variance of the prior on the constant
mn.psi = 1;   % =1 => maximizing also with respect to the scale parameters of the MN prior variance
mn.alpha = 0; % =1 => maximizing also with respect to the lag decaying parameter of the MN prior
fcast = 0;    % =1 => produces and evaluates the MSE of the forecasts
hz = 1;   % forecasts horizons

mcmc = 1;     % =1 => generates MCMC draws from the posterior
burn = 5000; % burning sample
MCMCfcast = 1;% =1 => generates density forecasts using the MCMC draws
hyperpriors = 1;    % hyperpriors = 0: no priors on hyperparameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load DATA_1970-2012_AWM_scale

Start = find(DateV(:,1)==StartE(1) & DateV(:,2)==StartE(2));
Last = find(DateV(:,1)==EndE(1) & DateV(:,2)==EndE(2));

% dates
DateV = DateV(Start:Last,:);
DateN = DateN(Start:Last,:);
DateS = DateS(Start:Last,:);

% differenced (log-)data
dy = dlogX(Start:Last,:);
% (log-)level data
y  = logX(Start:Last,:);

[T,nobs] = size(y);

% index in 1999 = 100
idx99 = find(DateV(:,1)==1999);
ScShift = mean(y(idx99,:),1)-log(100);
y = y - repmat(ScShift,T,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter draws
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% drawing the parameters using Giannone, Lenza and Primiceri procedure
yEst = y(1:end-nQoos,:); % estimation set, model in levels
dyEst = dy(1:end-nQoos,:); % estimation set, model in differences

% BVAR IN LEVELS
const = 0.5; %Constant in MCMC algorithm, should be calibrated to get 15 to 35%
           %acceptance rate
sur = 1;      % =1 => sur prior on (maximizing with respect to the variance)
noc = 1;      % =1 => noc prior on (maximizing with respect to the variance)
lags = p+1;     % # lags in the VAR
pos = [];     % position of variables whit prior centered on 0
seedL = rng('default');
%rL = optimalprior_mic(yEst,lags,Vc,pos,mn,sur,noc,fcast,hz,mcmc,nDraws+burn,MCMCfcast,hyperpriors,const);
rL2 = bvarGLP(dyEst, lags, 'MNpsi', mn.psi, 'MNalpha', mn.alpha, 'sur', sur,...
    'noc', noc, 'mcmc', mcmc, 'MCMCfcast', MCMCfcast, 'Ndraws', nDraws, ...
    'Ndrawsdiscard', burn, 'fcast', fcast, 'hz', hz, ...
    'hyperpriors', hyperpriors, 'MCMCconst', const, 'Vc', Vc);


% BVAR IN DIFFERENCES
const = 0.5; %Constant in MCMC algorithm, should be calibrated to get 15 to 35%
           %acceptance rate
sur = 0;      % =1 => sur prior on (maximizing with respect to the variance)
noc = 0;      % =1 => noc prior on (maximizing with respect to the variance)
lags = p;     % # lags in the VAR
pos = (1:nobs);     % position of variables with prior centered on 0
seedD = rng;

% rD = optimalprior_mic(dyEst, lags, Vc, pos, mn, sur, noc, fcast, hz, mcmc,...
%                       nDraws+burn, MCMCfcast, hyperpriors, const);    

rD2 = bvarGLP(dyEst, lags, 'MNpsi', mn.psi, 'MNalpha', mn.alpha, 'sur', sur,...
    'noc', noc, 'mcmc', mcmc, 'MCMCfcast', MCMCfcast, 'Ndraws', nDraws, ...
    'Ndrawsdiscard', burn, 'fcast', fcast, 'hz', hz, ...
    'hyperpriors', hyperpriors, 'MCMCconst', const, 'Vc', Vc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conditional forecasts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conditioning variables
idxCV = ismember(NewNames,CV);

% estimation and conditioning sets
dyT = dy;
dy(nQinit+1:end,~idxCV) = nan; % conditioning set, model in differences
yT = y;
y(nQinit+1:end,~idxCV) = nan; % conditioning set, model in levels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DFM

dydfm3 = DFMcf(dy,r,q,p,500,dyEst);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BVAR IN DIFFERENCES

% density conditional forecasts
dybvarD = nan(T,nobs, nDraws);
for i = 1:nDraws-burn
    clc
    disp(['Processing ' num2str(i) ' of ' num2str(nDraws - burn) '...'])
    tmp   = squeeze(rD2.mcmc.beta(:,:,i));
    Gamma = [tmp(2:end,:);tmp(1,:)];
    Su    = squeeze(rD2.mcmc.sigma(:,:,i));
    tmp   = VARcf_DKcks(dy,p,Gamma,Su,1);
    dybvarD(:,:,i) = tmp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BVAR IN LEVELS

% density conditional forecasts
ybvarL = nan(T,nobs, nDraws);
for i=1:nDraws-burn
    clc
    disp(['Processing ' num2str(i) ' of ' num2str(nDraws - burn) '...'])

    
    tmp   = squeeze(rL2.mcmc.beta(:,:,i));
    Gamma = [tmp(2:end,:);tmp(1,:)];
    Su    = squeeze(rL2.mcmc.sigma(:,:,i));
    tmp   = VARcf_DKcks(y, p+1, Gamma, Su, 1);
    ybvarL(:,:,i) = tmp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cumulate the results for the models in differences to get the levels
ybvarD        = dybvarD;
ybvarD(1,:,:) = repmat(yT(1,:),[1 1 nDraws]);
ybvarD        = cumsum(ybvarD,1);

ydfm3      = dydfm3;
ydfm3(1,:) = yT(1,:);
ydfm3      = cumsum(ydfm3,1);

% Variables to be plotted
ChoicePlot = {'YWR';'ITR';'XTR';'MTR';'URX';'PPIXC'; 'POILU';'LTN';'M1';'M3';'LOAH';'LOAC'};
% variables to be plotted in levels
LevelVars = {'URX';'LTN'};

% fan chart
% density conditional forecast from the level BVAR, median conditional forecast from the difference BVAR
% and point conditional forecast from the DFM
xaxislab = DateV(nQinit+1:end,1)+DateV(nQinit+1:end,2)/12-.25;
for s=1:12
    ser_s = ChoicePlot(s);
    idx_s = find(ismember(NewNames,ser_s));
    NewName_s = LongNames(idx_s);
%     figure
    subplot(3,4,s)
    tmp1 = squeeze(ybvarL(:,idx_s,:));
    tmp2 = squeeze(ybvarD(:,idx_s,:));
    tmp3 = [ydfm3(:,idx_s),yT(:,idx_s)];
    if ~ismember(ser_s,LevelVars)
        % annual growth rates
        tmp1 = cat(1,nan(4,nDraws),100*(tmp1(5:end,:)-tmp1(1:end-4,:))/4);
        tmp2 = cat(1,nan(4,nDraws),100*(tmp2(5:end,:)-tmp2(1:end-4,:))/4);
        tmp3 = cat(1,nan(4,2),100*(tmp3(5:end,:)-tmp3(1:end-4,:))/4);
    else
        tmp1 = 100*(tmp1+ScShift(idx_s));
        tmp2 = 100*(tmp2+ScShift(idx_s));
        tmp3 = 100*(tmp3+ScShift(idx_s));
    end
    tmp1 = sort(tmp1,2);
    tmp2 = sort(tmp2,2);    

    FanChart4(tmp1(nQinit+1:end,:),tmp2(nQinit+1:end,round(nDraws/2)),tmp3(nQinit+1:end,1),tmp3(nQinit+1:end,2),NewName_s,xaxislab)
end

save Res_CondFore
