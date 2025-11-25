%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figure 7

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed

saveFolder = 'Figures';
%% Load data

makedata_ea; % organizes the raw data into a T x n matrix and stores dates

%% SVAR with sign restrictions - SUR prior

draws = 10000; % Number of final draws
ndDisc = 10000; % Number of draws to discard

p = 4; % Lags
iid = [1,2]; % Set all variables to have 0 mean on own first lag
y0_bar = mean(data,1); % Chose y0_bar
[T,n] = size(data);

% Estimate reduced form model
res_ea = bvarGLP_y0(data, p, 'noc',0,'sur', 1, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws*1+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]

% Get posterior draws for beta and sigma
Bdraw_all        = res_ea.mcmc.beta;
Sigmadraw_all     = res_ea.mcmc.sigma;

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
plotDates = datetime(dates(p+1:end,:));


%% Fig 7 a: HD top 3
% Historical decompositions short - top 3

hFig = figure;
scale = 3;
numColumns = 3;
numRows = 1.7;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
titles = {'Closest','2nd closest','3rd closest'};

n_quarters = 12*3;
for k=1:3
    subplot(1,3,k)
    bar(plotDates(end-n_quarters:end),[squeeze(initialcond(end-n_quarters:end,2,solutions(k))) squeeze(histdec(end-n_quarters:end,:,2,solutions(k)))],'stacked'), hold on;
    plot(plotDates(end-n_quarters:end),Y(end-n_quarters:end,2),'LineWidth',2,'Color','k'), axis tight;
    
    title(titles{k})
    
    ylim([-4 12])
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

filename = strcat(saveFolder,'/Fig7a_fred');
print(hFig, filename, '-dpdf','-r0')

%% Fig 7b, left: All DC

tooHigh = find(initialcond(end,2,:)>6);
tooLow = find(initialcond(end,2,:)<-2);
initialcond_noExtremes = initialcond;
initialcond_noExtremes(:,:,[tooLow;tooHigh]) = [];
size(initialcond_noExtremes)
size(plotDates)
hFig = figure;
scale = 5;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])

i = 2;
plot(plotDates,squeeze(initialcond_noExtremes(:,i,:))), hold on;
plot(plotDates,Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;
ylim([-4 12])
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/Fig7bL_fred');
print(hFig, filename, '-dpdf','-r0')

%% Fig 7b. right: DC top 10

hFig = figure;
scale = 5;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
VARnames = {'GDP growth','Inflation','Fed funds rate'};
i = 2;

plot(plotDates,squeeze(initialcond(:,i,solutions))), hold on;
plot(plotDates,Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;
ylim([-4 12])
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 


set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/Fig7bR_fred');
print(hFig, filename, '-dpdf','-r0')



