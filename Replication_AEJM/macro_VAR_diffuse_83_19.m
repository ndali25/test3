%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figure A-3

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed
saveFolder = 'Figures';
priorName = 'diffuse_83_19';

%% Load data

makedata_83_19;

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

%% Figure A-3
VARnames={'GDP growth'; 'Inflation'};

tooHigh = find(initialcond(end,2,:)>2);
tooLow = find(initialcond(end,2,:)<-0.5);
initialcond_noExtremes = initialcond;
initialcond_noExtremes(:,:,[tooLow;tooHigh]) = [];

hFig = figure;
scale = 5;
numColumns = 2;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
for i = 1:2
    subplot(1,2,i)
    plot(plotDates,squeeze(initialcond_noExtremes(:,i,:))), hold on;
    plot(plotDates,Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;
    if i==1
        ylim([-2.5, 2.5])
    else
        ylim([-0.5, 2.5])
    end
    title(VARnames{i})
    set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
end
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/FigA3');
print(hFig, filename, '-dpdf','-r0')

