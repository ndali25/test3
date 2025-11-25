%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figure C-1

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed
saveFolder = 'Figures';

%% Simulate data for a VAR(1) with three variables:


% T=500; % Length of the data
T=150; % Length of the data
% T= 80; % Length of the data

n=2; % variables 

%A=[0.6 -0.3;0.3 0.4]; % Autoregressive matrix A
% A=[0.95 -0.3;0.3 0.4]; % Autoregressive matrix A
A=[0.2 -0.3;0.3 0.4]; % Autoregressive matrix
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

%% Estimate a VAR(1) with three variables using simulated data:

data = x;

p=1;
c=1;

% Draws set up:

draws=10000;

% Set up the loop for each draw :

PI=zeros(n*p+c,n,draws);
BigA=zeros(n*p,n*p,draws);
Sigma=zeros(n,n,draws);
errornorm=zeros(T-p,n,draws);
fittednorm=zeros(T-p,n,draws);

% Reduced-form VAR estimation - flat prior:
options.priors.name = 'Jeffrey';
options.K = draws;
BVAR = bvar_(data,p,options);
Bdraw_all        = [BVAR.Phi_draws(end,:,:); BVAR.Phi_draws(1:end-1,:,:)]; % Organize beta matrix with constant first
Sigmadraw_all     = BVAR.Sigma_draws;

% % Reduced-form VAR estimation - single-unit-root prior:
% ndDisc = 1000;
% mcmc_const = 4;
% y0_bar = mean(data,1);
% iid = [1,2];
% res    = bvarGLP_y0(data, p, 'noc',0,'sur', 1, 'mcmc', 1, 'MCMCconst', mcmc_const, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
% Bdraw_all = res.mcmc.beta;
% Sigmadraw_all = res.mcmc.sigma;

PI = Bdraw_all;
Sigma = Sigmadraw_all;
for i=1:draws 
    BigA(:,:,i)         = [PI(1+c:end,:,i)'; eye(n*p-n) zeros(n*p-n,n)]; % (n*p)x(n*p) matrix
%     Keep only stable draws
    if abs(eig(BigA(:,:,i)))>=1
        disp('Unstable')
        PI(:,:,i) = [];
        BigA(:,:,i)= [];
    end 
    [Y,X,Y_initial]     = SUR(data,p,c);
    errornorm(:,:,i)    = Y-X*PI(:,:,i);
    fittednorm(:,:,i)   = X*PI(:,:,i);
end

%% Historical decomposition and deterministic component:

% Data:

Y=x(p+1:end,:); % Discard first p observations

% Stochastic component:

err=errornorm; % Reduced-form residuals
HH=zeros(T-p,n,n,draws); % Store the impact IRFs ^t 
At=zeros(T-p,n,n,draws); % Store A^t 
hist=zeros(T-p,n,n,draws);

for k=1:draws
for t=1:T-p  
    
        BigC=(BigA(:,:,k))^(t-1);
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

for k=1:draws % Draws
for t=1:T-p % Time

    for i=1:n % Variables 
        
        initialcond(t,i,k)=Y(t,i)-sum(hist(t,:,i,k),2);  
        
    end   

    Aty0(t,:,k)=squeeze(At(t,:,:,k))*x(1,:)';

end
end

%% Plot initial conditions together with data:

VARnames={'$Y_1$';'$x_2$'};

hFig = figure;
for i=1:1
% subplot(n,1,i)
plot(squeeze(initialcond(:,i,:))), hold on;
plot(Y(:,i),'LineWidth',2,'Color','k'), axis tight;
title(VARnames{i},'Interpreter','Latex')
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
end

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = 'Figures/FigC1';
print(hFig, filename, '-dpdf','-r0')

