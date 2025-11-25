%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figure 3.(e)

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed

saveFolder = 'Figures';
priorName = 'diffuse_long_sample';

%% Load data

makedata_long; % organizes the raw data into a T x n matrix and stores dates

%% SVAR with sign restrictions - diffuse prior

draws = 10000; % Number of final draws
ndDisc = 0; % Number of draws to discard

p = 4; % Lags
c = 1; % constant

% Set up the loop for each draw :
PI=zeros(n*p+c,n,draws);
BigA=zeros(n*p,n*p,draws);
Sigma=zeros(n,n,draws);
errornorm=zeros(T-p,n,draws);
fittednorm=zeros(T-p,n,draws);
    
% Reduced-form VAR estimation:
    
for i=1:draws    
    stable=-1;    
    while stable<0
        
    [PI(:,:,i),BigA(:,:,i),Sigma(:,:,i),errornorm(:,:,i),fittednorm(:,:,i)]=BVAR(data,p,c);
        
    if abs(eig(BigA(:,:,i)))<1
            stable=1; % keep only stable draws
    end  
    end   
end

Bdraw_all = PI;
Sigmadraw_all = Sigma;

% Identify shocks
hor = 17; % IRF length in the plot

% Sign restricions
[f,sr] = var_restrictions('sign');
noConstant = 0;
identify_shocks;

% Compute historical decompositions
hist_decomp;

% Obtain the draws closest to the point-wise median IRF
fry_pagan;

%% Figures
plotDates = dates(p+1:end,:);

%% Figure 3.(e), left panel: All DC

tooHigh = find(initialcond(end,2,:)>2);
tooLow = find(initialcond(end,2,:)<-0.5);
initialcond_noExtremes = initialcond;
initialcond_noExtremes(:,:,[tooLow;tooHigh]) = [];

hFig = figure;
scale = 5;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])

i = 2;
plot(plotDates,squeeze(initialcond_noExtremes(:,i,:))), hold on;
plot(plotDates,Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;

set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off');

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/Fig3e');
print(hFig, filename, '-dpdf','-r0')

