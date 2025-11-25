%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figure 8 (a) and (b) - Australia 

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed

saveFolder = 'Figures';
%% Load data
datafile = 'Data/IMF_data_different_countries.xlsx';
[num,txt] = xlsread(datafile,'AU');
dates_raw = txt(2:end,1);
startDate = '01.01.1993';
endDate = '01.04.2023';
start_ind = find(strcmp(startDate,dates_raw));
end_ind = find(strcmp(endDate,dates_raw));
data_excel = num(start_ind:end_ind,:);
raw_variables = txt(1,2:end);
dates = datetime(dates_raw(start_ind:end_ind));

for ii = 1:length(raw_variables)
    rawData.(raw_variables{ii}) = data_excel(:,ii);
end

%% Convert data

variable_list = {%  Short name, code generating variable,           name,           transformation
                    'y',        'log(rawData.NGDP_R_SA_XDC)',       'GDP',          'diff';...
                    'p',        'log(rawData.PCPI_IX)',             'Inflation',    'yoy'};


variables = variable_list(:,1)'; % Short name of variable
var_legends = variable_list(:,3)'; % Descriptive name of variable

% Load data into the DATA object
for ii = 1:size(variable_list,1)
    codeString = strcat('DATA.',variable_list{ii,1},' = ',variable_list{ii,2},';');
    eval(codeString);
end

data    = nan(size(DATA.y,1)-4,size(variables,2)); % Data matrix
DATA
% Transform data
for jj = 1:size(variables,2)
    if strcmp('diff',variable_list{jj,4})
        data(:,jj) = DATA.(variables{jj})(5:end,:)-DATA.(variables{jj})(4:end-1,:);
    elseif strcmp('yoy',variable_list{jj,4})
        data(:,jj) = DATA.(variables{jj})(5:end,:)-DATA.(variables{jj})(1:end-4,:);
    else
        data(:,jj) = DATA.(variables{jj})(5:end,:);
    end
end

% Multiply by 100
data = data*100;

dates = dates(2:end);

[T,n] = size(data);
%% General settings
draws = 10000; % Number of final draws
ndDisc = 5000; % Number of draws to discard


saveFolder = 'Figures';
priorName = 'sur_AU';
%% SVAR with sign restrictions - SUR prior

p = 4; % Lags
iid = [1,2]; % Set all variables to have 0 mean on own first lag
y0_bar = mean(data,1); % Chose y0_bar

% Estimate reduced form model
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
plotDates = dates(p+4:end,:);


%% Fig 6 a 4: HD top 3
% Historical decompositions short - top 3

hFig = figure;
scale = 3;
numColumns = 1;
numRows = 1.7;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
titles = {'Australia'};

n_quarters = 28;
for k=1:1
% subplot(1,3,k)
bar(plotDates(end-n_quarters:end),[squeeze(initialcond(end-n_quarters:end,2,solutions(k))) squeeze(histdec(end-n_quarters:end,:,2,solutions(k)))],'stacked'), hold on;
plot(plotDates(end-n_quarters:end),Y(end-n_quarters:end,2),'LineWidth',2,'Color','k'), axis tight;

title(titles{k})

ylim([-6, 8])
pos = get(gca,'Position');
pos(2) = pos(2)+.1;
pos(4) = pos(4)-.2;
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
end
% legend({'Deterministic component','Supply','Demand'},'Orientation','Horizontal','Position',[0.5 0.02 0 0],'Box','off')

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/Fig8a4_HD_FP1_',priorName);
print(hFig, filename, '-dpdf','-r0')



%% Fig 8b 4: DC top 10

hFig = figure;
scale = 5;
numColumns = 1;
numRows = 1;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
VARnames = {'GDP growth','Inflation'};
i = 2;

plot(plotDates,squeeze(initialcond(:,i,solutions))), hold on;
plot(plotDates,Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;
ylim([-6, 8])
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 


set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/Fig8b4_DC_all_',priorName);
print(hFig, filename, '-dpdf','-r0')






