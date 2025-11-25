%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figure D-1.(c)

clear; clc;

MTBfolder = 'Functions/';

addpath(strcat(MTBfolder,'getfreddata-matlab-master'))
addpath(strcat(MTBfolder,'Functions_Primiceri/GLPreplicationWeb'))
addpath(strcat(MTBfolder,'Functions_Primiceri/BGLreplication'))
addpath(strcat(MTBfolder,'Functions_Primiceri/GLPreplicationWeb/subroutines'))
addpath(strcat(MTBfolder,'canova/BVAR_-master/bvartools'))
addpath(strcat(MTBfolder,'matlab_functions'))
addpath(strcat(MTBfolder,'binning'))

set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

seed=12345;
rng(seed)

%% Simulate data for a VAR(1) with three variables:

T=10000; % Length of the data

n=2; % variables 

A=[1 0;0 1]; % Autoregressive matrix A

C=[0;0];
%B0=[0.6 0;-0.3 0.5]; % Structural matrix B0
B0=eye(2); % Structural matrix B0

e=standardize(randn(T,n)); % Generate structural errors N(0,1)

x(:,1)=zeros(n,1)+randn(n,1); % Data
u(:,1)=zeros(n,1); % Reduced-form errors

for t=2:T
        
        u(:,t)=B0*e(t,:)';
        x(:,t)=C+A*x(:,t-1)+u(:,t);
        
end

x = x';
u = u';

%% Estimate a VAR(1) with three variables using simulated data:

code = 1;
vardata = x;

minnesota = 1;
sur = 1;
noc = 0;

lambda = []; %0.2;
theta  = []; % 0.0001;
mu     = [];%1;

iid = [1,2];

p=1;
c=1;

y0_bar = mean(vardata,1);
% Draws set up:

draws=1000;

% Set up the loop for each draw :

PI=zeros(n*p+c,n,draws);
BigA=zeros(n*p,n*p,draws);
Sigma=zeros(n,n,draws);
errornorm=zeros(T-p,n,draws);
fittednorm=zeros(T-p,n,draws);

% Reduced-form VAR estimation - flat prior:

if code==1
    % Reduced-form VAR estimation:

    for i=1:draws    
        stable=-1;    
        while stable<0

        [PI(:,:,i),BigA(:,:,i),Sigma(:,:,i),errornorm(:,:,i),fittednorm(:,:,i)]=BVAR(vardata,p,c);

        if abs(eig(BigA(:,:,i)))<1
            stable=1; % keep only stable draws
        end  
        end   
    end
elseif code==2
    ndDisc = 1000;
    mcmc_const = 4;
    res    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 1, 'MCMCconst', mcmc_const, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
    PI = res.mcmc.beta;
    Sigma = res.mcmc.sigma;
    for i=1:draws 
        [Y,X,Y_initial]     = SUR(vardata,p,c);
        BigA(:,:,i)         = [PI(1+c:end,:,i)'; eye(n*p-n) zeros(n*p-n,n)]; % (n*p)x(n*p) matrix
        if abs(eig(BigA(:,:,i)))>=1
            disp('Unstable')
            PI(:,:,i) = [];
            BigA(:,:,i)= [];
        end  
        errornorm(:,:,i)    = Y-X*PI(:,:,i);
        fittednorm(:,:,i)   = X*PI(:,:,i);
    end
    
elseif code==5
    
    y0_bar = mean(vardata,1);
    iid = [1,2];
    res    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 0, 'MNalpha',0,'MNpsi',0,'pos', iid, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
    

    lambda = []; %0.2;
    theta  = [0.0001]; % 0.0001;
    mu     = [];%1;
    
    prior_settings.prior = 'Minnesota';

    prior_settings.Minnesota = minnesota;
    prior_settings.SingleUnitRoot = sur;
    prior_settings.NoCointegration = noc;

    prior_settings.onlyprior = 0;

    prior_settings.stationary = [];
    prior_settings.iid = iid;
    if minnesota==1 && isempty(lambda)
        prior_settings.hyperparameters.lambda = res.postmax.lambda;
    elseif minnesota==1 && ~isempty(lambda)
        prior_settings.hyperparameters.lambda = lambda;
    end
    if sur==1 && isempty(theta)
        prior_settings.hyperparameters.theta = res.postmax.theta;
    elseif sur==1 && ~isempty(theta)
        prior_settings.hyperparameters.theta = theta;
    end
    if noc==1 && isempty(mu)
        prior_settings.hyperparameters.mu = res.postmax.miu;
    elseif noc==1 && ~isempty(mu)
        prior_settings.hyperparameters.mu = mu;
    end

    constant = c;       % Include Intercept
    y = vardata;
    n=size(y,2);

    prior_settings.hyperparameters.Y0bar = y0_bar;

    VAR_rubio
    Bdraw_all = B_draws;
    Sigmadraw_all = SIGMA_draws;
    

T=150;    
Bdraw_all_cFirst = Bdraw_all;

BigA=zeros(n*p,n*p,draws);
Sigma=zeros(n,n,draws);
errornorm=zeros(T-p,n,draws);
fittednorm=zeros(T-p,n,draws);
for i=1:draws 
    [Y,X,Y_initial]     = SUR(vardata,p,c);
    BigA(:,:,i)         = [Bdraw_all_cFirst(1+c:end,:,i)'; eye(n*p-n) zeros(n*p-n,n)]; % (n*p)x(n*p) matrix
    if abs(eig(BigA(:,:,i)))>=1
        disp('Unstable')
        Bdraw_all_cFirst(:,:,i) = [];
        BigA(:,:,i)= [];
    end  
    errornorm(:,:,i)    = Y-X*Bdraw_all_cFirst(:,:,i);
    fittednorm(:,:,i)   = X*Bdraw_all_cFirst(:,:,i);
end
    
end

%% Historical decomposition and initial conditions:

% Data:

Y=x(p+1:end,:); % Discard first p observations

% Stochastic component:

err=errornorm; % Reduced-form residuals
HH=zeros(T-p,n,n,draws); % Store the impact IRFs ^t 
At=zeros(T-p,n,n,draws); % Store A^t 
hist=zeros(T-p,n,n,draws);
InvA=zeros(n*p,n*p,draws);

for k=1:draws
for t=1:T-p  
    
        BigC=(BigA(:,:,k))^(t-1);
        InvA(:,:,k)=inv(eye(n*p,n*p)-BigA(:,:,k));
        BigH=BigC(1:n,1:n);
        At(t,:,:,k)=BigH;
        
        for j=1:n
            
        HH(t,:,j,k)=BigH(j,1:n);
        
        end       
end

%HH=cumsum(HH,1); % If variables are in first differences

for t=1:T-p
    for i=1:n
        for j=1:n
        
          hist(t,i,j,k)=HH(1:t,i,j,k)'*flipud(err(1:t,i,k));
          
        end        
    end   
end
end

% Deterministic component:

initialcond=zeros(T-p,n,draws);
Aty0=zeros(T-p,n,draws);
Aty0dev=zeros(T-p,n,draws);
sumAC=zeros(T-p,n,draws);
lrmean=zeros(T-p,n,draws);

for k=1:draws % Draws
for t=1:T-p % Time

    for i=1:n % Variables 
        
        initialcond(t,i,k)=Y(t,i)-sum(hist(t,:,i,k),2);  
        
    end   

    Aty0(t,:,k)=squeeze(At(t,:,:,k))*x(1,:)';
    Aty0dev(t,:,k)=squeeze(At(t,:,:,k))*(x(1,:)'-InvA(1:n,1:n,k)*PI(1,:,k)');
    sumAC(t,:,k)=squeeze(At(t,:,:,k))*PI(1,:,k)';
    lrmean(t,:,k)=InvA(1:n,1:n,k)*PI(1,:,k)';


end
end

%% Plot initial conditions together with data:

tooHigh = find(initialcond(end,1,:)>prctile(initialcond(end,1,:),95));

tooLow = find(initialcond(end,1,:)<prctile(initialcond(end,1,:),5));
initialcond_noExtremes = initialcond;
initialcond_noExtremes(:,:,[tooLow;tooHigh]) = [];
Aty0dev_noExtremes=Aty0dev;
Aty0dev_noExtremes(:,:,[tooLow;tooHigh]) = [];
lrmean_noExtremes=lrmean;
lrmean_noExtremes(:,:,[tooLow;tooHigh]) = [];



hFig = figure;
scale = 5;
numColumns = 3;
numRows = 2;

set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])

subplot(1,3,1)
plot(squeeze(initialcond_noExtremes(:,1,:))), hold on;
plot(Y(:,1),'LineWidth',2,'Color','k','LineStyle','-'), axis tight;
title('Deterministic component','Interpreter','Latex')
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

subplot(1,3,2)
plot(squeeze(lrmean_noExtremes(:,1,:))), hold on;
plot(Y(:,1),'LineWidth',2,'Color','k','LineStyle','-'), axis tight;
title('$(\mathbf{I}-\mathbf{A})^{-1}\mathbf{C}$','Interpreter','Latex')
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

subplot(1,3,3)
plot(squeeze(Aty0dev_noExtremes(:,1,:))), hold on;
plot(Y(:,1),'LineWidth',2,'Color','k','LineStyle','-'), axis tight;
title('$\mathbf{A}^t(\mathbf{Y}_0-(\mathbf{I}-\mathbf{A})^{-1}\mathbf{C})$','Interpreter','Latex')
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

filename = 'Figures/FigD1c';
print(hFig, filename, '-dpng','-r0')