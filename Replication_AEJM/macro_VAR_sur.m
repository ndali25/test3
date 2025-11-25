%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figures 6, 10.(d) and A-2

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed

saveFolder = 'Figures';
priorName = 'sur';

%% Load data

makedata; % organizes the raw data into a T x n matrix and stores dates

%% SVAR with sign restrictions - SUR prior with sample mean as y0_bar

draws = 10000; % Number of final draws
ndDisc = 3000; % Number of draws to discard

p = 4; % Lags
iid = [1,2]; % Set all variables to have 0 mean on own first lag
y0_bar = mean(data,1); % Chose y0_bar

% Estimate reduced form model with sur prior
res_sur = bvarGLP_y0(data, p, 'noc',0,'sur', 1, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws*1+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]

% Get posterior draws for beta and sigma
Bdraw_all        = res_sur.mcmc.beta;
Sigmadraw_all     = res_sur.mcmc.sigma;

% Identify shocks
hor = 17; % IRF length
% Set restricions
[f,sr] = var_restrictions('sign');
noConstant = 0;
identify_shocks;

% Compute historical decompositions
hist_decomp;

% Obtain the draws closest to the point-wise median IRF
fry_pagan;

%% Figures
plotDates = dates(p+1:end,:);

%% Fig 6.(a): HD top 3
% Historical decompositions short - top 3

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
plot(plotDates(end-n_quarters:end),Y(end-n_quarters:end,2),'LineWidth',2,'Color','k'), axis tight;

title(titles{k})

ylim([-2, 2.5])
pos = get(gca,'Position');
pos(2) = pos(2)+.1;
pos(4) = pos(4)-.2;
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
end
legend({'Deterministic component','Supply','Demand'},'Orientation','Horizontal','Position',[0.5 0.02 0 0],'Box','off')

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/Fig6a');
print(hFig, filename, '-dpdf','-r0')

%% Fig 6.(b), left panel: All DC

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
ylim([-0.5, 2.5])
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/Fig6bL');
print(hFig, filename, '-dpdf','-r0')

%% Fig 6.(b), right panel: DC top 10

hFig = figure;
scale = 5;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
VARnames = {'GDP growth','Inflation','Fed funds rate'};
i = 2;

plot(plotDates,squeeze(initialcond(:,i,solutions))), hold on;
plot(plotDates,Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;
ylim([-0.5, 2.5])
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 


set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/Fig6bR');
print(hFig, filename, '-dpdf','-r0')

%% Fig 10.(d): Median historical decomposition - 2019 onwards

hFig = figure;
scale = 5;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])

n_quarters = 15;
residual_component = Y(end-n_quarters:end,2)-sum(squeeze(median(histdec(end-n_quarters:end,:,2,:),4)),2);


bar(plotDates(end-n_quarters:end),[residual_component squeeze(median(histdec(end-n_quarters:end,:,2,:),4))],'stacked'), hold on;
plot(plotDates(end-n_quarters:end),Y(end-n_quarters:end,2),'LineWidth',2,'Color','k');
legend({'Residual component','Supply','Demand'},'Orientation','Horizontal','Location', 'southoutside','Box','off') 
ylim([-2, 2.5])
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
filename = strcat(saveFolder,'/Fig10d');
print(hFig, filename, '-dpdf','-r0')


%% Fig A-2: Delta
hFig = figure;

% Plot the histogram
histogram(res_sur.mcmc.theta, 10,"BinLimits",[0 0.0005]); % Or use binEdges instead of numBins
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off');

% hist(res_ea.mcmc.theta)
title('Posterior of the shrinckage of the single-unit-root prior');
mean(res_sur.mcmc.theta)

filename = strcat(saveFolder,'/FigA2');
print(hFig, filename, '-dpdf','-r0')
