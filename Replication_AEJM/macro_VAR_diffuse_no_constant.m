%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figure 9

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(12345) % Set seed

saveFolder = 'Figures';
priorName = 'diffuse_no_constant';

%% Load data

makedata; % organizes the raw data into a T x n matrix and stores dates

%% SVAR with sign restrictions - diffuse prior with de-meaned data

draws = 10000; % Number of final draws
ndDisc = 0; % Number of draws to discard

% De-mean data
data = data-mean(data);

p = 4; % Lags
c = 0; % constant

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
noConstant = 1;
identify_shocks;

% Compute historical decompositions
hist_decomp;

% Obtain the draws closest to the point-wise median IRF
fry_pagan;

%% Figures
plotDates = dates(p+1:end,:);

%% Figure 9.(a): HD top 3

% Historical decompositions short - top 3

yLimit = [-2, 2.5];
hFig = figure;
scale = 3;
numColumns = 3;
numRows = 1.7;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
titles = {'Closest','2nd closest','3rd closest'};

n_quarters = 28;
for k=1:3
subplot(1,3,k)
bar(plotDates(end-n_quarters:end),[squeeze(initialcond(end-n_quarters:end,2,solutions(k))) squeeze(histdec(end-n_quarters:end,:,2,solutions(k)))],'stacked'), hold on;
plot(plotDates(end-n_quarters:end),Y(end-n_quarters:end,2),'LineWidth',2.5,'Color','k'), axis tight;
title(titles{k})

ylim(yLimit)
pos = get(gca,'Position');
pos(2) = pos(2)+.1;
pos(4) = pos(4)-.2;
set(gca,'FontSize',18,'Position',pos)
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
end
legend({'Deterministic component','Supply','Demand'},'Orientation','Horizontal','Position',[0.5 0.05 0 0],'Box','off')

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% Save figure

filename = strcat(saveFolder,'/Fig9a');
print(hFig, filename, '-dpdf','-r0')

%% Figure 9.(b), left panel: All DC

tooHigh = find(initialcond(end,2,:)>2); 
tooLow = find(initialcond(end,2,:)<-0.5);
initialcond_noExtremes = initialcond;
initialcond_noExtremes(:,:,[tooLow;tooHigh]) = [];  % remove extreme DC

hFig = figure;
scale = 5;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])

i = 2;
plot(plotDates,squeeze(initialcond_noExtremes(:,i,:))), hold on;
plot(plotDates,Y(:,i),'-','LineWidth',2.5,'Color','k'), axis tight;
title('All posterior draws')

set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% Save figure

filename = strcat(saveFolder,'/Fig9bL');
print(hFig, filename, '-dpdf','-r0')

%% Figure 10.(b), right panel: DC top 10

hFig = figure;
scale = 5;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
VARnames = {'GDP growth','Inflation','Fed funds rate'};
i = 2;

plot(plotDates,squeeze(initialcond(:,i,solutions))), hold on;
plot(plotDates,Y(:,i),'-','LineWidth',2.5,'Color','k'), axis tight;
title('Top 10 draws')
ylim([-0.5, 2.5])
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 


set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% Save figure

filename = strcat(saveFolder,'/Fig9bR');
print(hFig, filename, '-dpdf','-r0')
