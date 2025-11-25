%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figure

clear; close all; clc;

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

%% Simulate data for a VAR(1) with two variables:
clear;
rng(12345)

T=150; % Length of the data (the 100 first observations are used as pre-sample)


n=2; % variables 

A=[0.95 -0.3;0.3 0.4]; % Autoregressive matrix A
C=[0.4;0.5];
B0=[0.6 0;-0.3 0.5]; % Structural matrix B0

e=standardize(randn(T,n)); % Generate structural errors N(0,1)

x(:,1)=zeros(n,1)+C; % Data
u(:,1)=zeros(n,1); % Reduced-form errors

for t=2:T
        
        u(:,t)=e(t,:)';
        x(:,t)=C+A*x(:,t-1)+u(:,t);
        
end

x = x';
u = u';

uncond_mean = (eye(n)-A)^(-1)*C;


%% Estimate a VAR(1) with two variables using simulated data:

code = 1;

vardata = x;

iid = [1,2];

p=1;
c=1;

delta_loose = 0.2;
delta_tight = 0.0001;

delta_all = [delta_loose, delta_tight];

disp_table = nan(size(delta_all,2),2,5);


y0_bar = mean(vardata,1);

T = size(vardata,1);
% Set up the loop for each draw :
T = T+1;
draws=1000;

PI=zeros(n*p+c,n,draws);
BigA=zeros(n*p,n*p,draws);
Sigma=zeros(n,n,draws);
errornorm=zeros(T-p,n,draws);
fittednorm=zeros(T-p,n,draws);

for d_i = 1:size(delta_all,2)
    delta = delta_all(d_i);
    for i=1:draws    
        stable=-1;    
        while stable<0

        [PI(:,:,i),BigA(:,:,i),Sigma(:,:,i),errornorm(:,:,i),fittednorm(:,:,i)]=BVAR_sur(vardata,p,c,y0_bar,delta);

        if abs(eig(BigA(:,:,i)))<1
            stable=1; % keep only stable draws
        end  
        end   
    end


    %% Historical decomposition and initial conditions:
    
    % Data:
    
    Y=x(p+1:end,:); % Discard first p observations
    
    % Stochastic component:
    T = T-1;
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
    
    if d_i == 1
        initialcond_loose = initialcond;
    elseif d_i ==2
        initialcond_tight = initialcond;
    end
   
end

%% Figure 5a: All DC, delta = 0.2

tooHigh = find(initialcond_loose(end,2,:)>2);
tooLow = find(initialcond_loose(end,2,:)<-0.5);
initialcond_noExtremes = initialcond_loose;
initialcond_noExtremes(:,:,[tooLow;tooHigh]) = [];

hFig = figure;
scale = 5;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])

i = 2;
plot(squeeze(initialcond_noExtremes(:,i,:))), hold on;
plot(Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;

set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = 'Figures/Fig5a';
print(hFig, filename, '-dpdf','-r0')

%% Figure 5b: All DC, delta = 0.0001

tooHigh = find(initialcond_tight(end,2,:)>2);
tooLow = find(initialcond_tight(end,2,:)<-0.5);
initialcond_noExtremes = initialcond_tight;
initialcond_noExtremes(:,:,[tooLow;tooHigh]) = [];

hFig = figure;
scale = 5;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])

i = 2;
plot(squeeze(initialcond_noExtremes(:,i,:))), hold on;
plot(Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;

set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = 'Figures/Fig5b';
print(hFig, filename, '-dpdf','-r0')

