%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figures 6, 10.(d) and A-2

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

if exist('seed','var') == 0
    seed = 7;
end
rng(seed)

%% Load data
datafile = 'Data/US_large.xlsx';
[num,txt] = xlsread(datafile);
dates_raw = txt(2:end,1);
startDate = '1983-01-01';
start_ind = find(strcmp(startDate,dates_raw));
data_excel = num(start_ind:end,:);
raw_variables = txt(1,2:end);
dates = datetime(dates_raw(start_ind:end));

for ii = 1:length(raw_variables)
    rawData.(raw_variables{ii}) = data_excel(:,ii);
end
%% Convert data

variable_list = {% Short name, code generating variable, name, transformation
    'y',    'log(rawData.GDPC1)',                                      'GDP',              'diff';...
    'p',    'log(rawData.GDPDEF)',                                    'Prices',           'diff';...
    'i',    'log(rawData.GPDIC1)',                                     'Investment',       'diff';...
    'r',    'rawData.FEDFUNDS*0.01',                                   'Fed funds rate',   'level';...
    'w',    'log(rawData.AHETPI)-log(rawData.GDPDEF)',              'Real wage',        'diff'};

variables = variable_list(:,1)'; % Short name of variable
var_legends = variable_list(:,3)'; % Descriptive name of variable

% Load data into the DATA object
for ii = 1:size(variable_list,1)
    codeString = strcat('DATA.',variable_list{ii,1},' = ',variable_list{ii,2},';');
    eval(codeString);
end

data    = nan(size(DATA.y,1)-1,size(variables,2)); % Data matrix
% dates = x2mdate(rawData.GDPC1.Data(2:end,1)-693960,0,'datetime');

% Take first differences (where specified)
for ii = 1:size(variables,2)
    if strcmp('diff',variable_list{ii,4})
        data(:,ii) = DATA.(variables{ii})(2:end,:)-DATA.(variables{ii})(1:end-1,:);
    else
        data(:,ii) = DATA.(variables{ii})(2:end,:);
    end
end

data = data*100;


dates = dates(2:end);
[T,n] = size(data);
%% General settings
draws = 1000; % Number of final draws
ndDisc = 1000; % Number of draws to discard


saveFolder = 'Figures';
priorName = 'sur_large_model';


%% SUR 83 with sample mean as y0_bar
rng(100) % Set seed
p = 4; % Lags
iid = [1,2,3,4,5]; % Set all variables to have 0 mean on own first lag
y0_bar = mean(data,1); % Chose y0_bar

% Estimate reduced form model
res_sur = bvarGLP_y0(data, p, 'noc',0,'sur', 1, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws*1+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]

% Get posterior draws for beta and sigma
Bdraw_all        = res_sur.mcmc.beta;
Sigmadraw_all     = res_sur.mcmc.sigma;

% Identify shocks
hor = 17; % IRF length
% Set restricions

noConstant = 0;
% Idenfify shocks

    %df, inv, mp, ls, other supply
f = [1,1,1,1,1; %GDP
     1,1,1,1,1; %Prices
     1,1,1,1,1; %Investment
     1,1,1,1,1; %Fed funds rate
     1,1,1,1,1; %Real wage
     %=========
     1,1,1,1,1; %GDP
     1,1,1,1,1; %Prices
     1,1,1,1,1; %Investment
     1,1,1,1,1; %Fed funds rate
     1,1,1,1,1]; %Real wage


      %df, inv, mp, ls, other supply
sr = [+1,+1,+1,+1,+1; %GDP
     +1,+1,+1,-1,-1; %Prices
     +1,+1,+1,nan,nan; %Investment
     +1,+1,-1,nan,nan; %Fed funds rate
     nan,nan,nan,-1,+1]; %Real wage

magnitude=1;

VARnames={'GDP'; 'Inflation';'Investment';'Fed funds rate';'Real wage'};
Shocknames={'Other demand';'Investment';'Monetary policy';'Labor supply';'Other supply'};

% draws = 100; % Number of draws
identify_shocks_large;

% Compute historical decompositions
hist_decomp;

% Obtain the draws closest to the point-wise median IRF
fry_pagan;


%% Figures
plotDates = datetime(dates(p+1:end,:));


%% Fig E-2: Median historical decomposition

% Historical decompositions
med_histDecInfl = median(squeeze(histdec(:,:,2,:)),3);
med_initialcond = median(squeeze(initialcond(:,2,:)),2);
residual_component = Y(:,2)-sum(med_histDecInfl,2);


hFig = figure;
scale = 6;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
% titles = {'Closest','2nd closest','3rd closest'};

n_quarters = 28;
% bar(plotDates(end-n_quarters:end),[med_initialcond(end-n_quarters:end,:) med_histDecInfl(end-n_quarters:end,:)],'stacked'), hold on;
bar(plotDates(end-n_quarters:end),[residual_component(end-n_quarters:end,:) med_histDecInfl(end-n_quarters:end,:)],'stacked'), hold on;
plot(plotDates(end-n_quarters:end),Y(end-n_quarters:end,2),'LineWidth',2,'Color','k'), axis tight;
ylim([-1.5 2.5])
% title(titles{k})
pos = get(gca,'Position');
pos(2) = pos(2)+.1;
pos(4) = pos(4)-.1;
set(gca,'FontSize',14,'Position',pos)
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

legend(['Residual component',Shocknames'],'Orientation','Horizontal','Position',[0.5 0.05 0 0],'Box','off')

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat('Figures/FigE2');
print(hFig, filename, '-dpdf','-r0')

%% Fig E-3: Mean historical decomposition

% Historical decompositions
mean_histDecInfl = mean(squeeze(histdec(:,:,2,:)),3);
mean_initialcond = mean(squeeze(initialcond(:,2,:)),2);
% residual_component = Y(:,2)-sum(med_histDecInfl,2);


hFig = figure;
scale = 6;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
% titles = {'Closest','2nd closest','3rd closest'};

n_quarters = 28;
bar(plotDates(end-n_quarters:end),[mean_initialcond(end-n_quarters:end,:) mean_histDecInfl(end-n_quarters:end,:)],'stacked'), hold on;
% bar(plotDates(end-n_quarters:end),[residual_component(end-n_quarters:end,:) med_histDecInfl(end-n_quarters:end,:)],'stacked'), hold on;
plot(plotDates(end-n_quarters:end),Y(end-n_quarters:end,2),'LineWidth',2,'Color','k'), axis tight;
ylim([-1.5 2.5])
% title(titles{k})
pos = get(gca,'Position');
pos(2) = pos(2)+.1;
pos(4) = pos(4)-.1;
set(gca,'FontSize',13,'Position',pos)
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

legend(['Deterministic component',Shocknames'],'Orientation','Horizontal','Position',[0.5 0.05 0 0],'Box','off')

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat('Figures/FigE3');
print(hFig, filename, '-dpdf','-r0')
