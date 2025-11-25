%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figures 3.(d) and 10.(c)

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed

saveFolder = 'Figures';
priorName = 'minnesota';

%% Load data

makedata; % organizes the raw data into a T x n matrix and stores dates

%% SVAR with sign restrictions - minnesota prior

draws = 10000; % Number of final draws
ndDisc = 3000; % Number of draws to discard

p = 4; % Lags
iid = [1,2]; % Set all variables to have 0 mean on own first lag
y0_bar = mean(data,1); % Chose y0_bar

% Estimate reduced form model
res_minnesota = bvarGLP_y0(data, p, 'noc',0,'sur', 0, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws*1+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]

% Get posterior draws for beta and sigma
Bdraw_all        = res_minnesota.mcmc.beta;
Sigmadraw_all     = res_minnesota.mcmc.sigma;

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

%% Fig 3.(d): All DC

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
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% Save figure

filename = strcat(saveFolder,'/Fig3d');
print(hFig, filename, '-dpdf','-r0')

%% Fig 10.(c): Median historical decomposition - 2019 onwards

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
filename = strcat(saveFolder,'/Fig10c');
print(hFig, filename, '-dpdf','-r0')

%% Uncomment here for IRFs and HD for minnesota prior

% %% Fig 6 a: HD top 3
% % Historical decompositions short - top 3
% 
% hFig = figure;
% scale = 3;
% numColumns = 3;
% numRows = 1.7;
% set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
% titles = {'Closest','2nd closest','3rd closest'};
% 
% n_quarters = 28;
% for k=1:3
% subplot(1,3,k)
% bar(plotDates(end-n_quarters:end),[squeeze(initialcond(end-n_quarters:end,2,solutions(k))) squeeze(histdec(end-n_quarters:end,:,2,solutions(k)))],'stacked'), hold on;
% plot(plotDates(end-n_quarters:end),Y(end-n_quarters:end,2),'LineWidth',2,'Color','k'), axis tight;
% 
% title(titles{k})
% 
% ylim([-2, 2.5])
% pos = get(gca,'Position');
% pos(2) = pos(2)+.1;
% pos(4) = pos(4)-.2;
% set(gca,'FontSize',18)
% end
% legend({'Deterministic component','Supply','Demand'},'Orientation','Horizontal','Position',[0.5 0.02 0 0],'Box','off')
% 
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% filename = strcat(saveFolder,'/FigA3b_HD_FP3_',priorName);
% print(hFig, filename, '-dpdf','-r0')
% 
% 
% %% IRF
% 
% conf=68;
% colorBNDS=[0.7 0.7 0.7];
% VARnames={'GDP growth'; 'Inflation'};
% Shocknames={'Aggregate supply';'Aggregate demand'};
% 
% 
% % 3 closest: Max/min
% hFig = figure;
% scale = 4;
% numColumns = n;
% numRows = n;
% set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
% 
% 
% 
% 
% 
% for k=1:n % variable k
%     for j=1:n % shock j
% 
%     subplot(n,n,j+n*k-n)
%     fill([0:hor-1 fliplr(0:hor-1)]' ,[reshape(permute(max(candidateirf(k,j,:,:),[],4),[3 2 1]),hor,1,[]); flipud(reshape(permute(min(candidateirf(k,j,:,:),[],4),[3 2 1]),hor,1,[]))],...
%         [.93 .93 .93],'EdgeColor','None'); hold on;
%     fill([0:hor-1 fliplr(0:hor-1)]' ,[reshape(permute(prctile(candidateirf(k,j,:,:),(100+conf)/2,4),[3 2 1]),hor,1,[]); flipud(reshape(permute(prctile(candidateirf(k,j,:,:),(100-conf)/2,4),[3 2 1]),hor,1,[]))],...
%         [.6 .6 .6],'EdgeColor','None'); hold on;
% %     plot(0:hor-1,reshape(permute(prctile(candidateirf(k,j,:,:),(100-conf)/2,4),[3 2 1]),hor,1,[]),'LineWidth',1.5,'Color','k'); hold on;
%     plot(0:hor-1,reshape(permute(prctile(candidateirf(k,j,:,:),50,4),[3 2 1]),hor,1,[]),'LineWidth',3.5,'Color','k'); hold on;
%     plot(0:hor-1,candidateirf_wold_fp(:,j+n*k-n),'LineWidth',3.5,'Color','r'); hold on;
%     plot(0:hor-1,candidateirf_wold_fp_2(:,j+n*k-n),'LineWidth',3.5,'Color','g'); hold on;
%     plot(0:hor-1,candidateirf_wold_fp_3(:,j+n*k-n),'LineWidth',3.5,'Color','b'); hold on;
% %     plot(0:hor-1,reshape(permute(prctile(candidateirf(k,j,:,:),(100+conf)/2,4),[3 2 1]),hor,1,[]),'LineWidth',1.5,'Color','k'); hold on;
%     line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5); hold off;
%     set(gca,'FontSize',16)
% 
%     if j ==1
%         ylabel(strcat(VARnames{k}), 'FontSize', 18, 'interpreter', 'latex');
%     end
% 
%     if k == 1
%         title(Shocknames{j}, 'FontSize', 16, 'interpreter', 'latex');
%     end
% 
%     xlim([0 hor-1]);
% 
%     end
% end
% legend({'min/max','68% credible set','Pointwise median','Fry-Pagan draw','2nd closest','3rd closest'},'FontSize',16,'Orientation','Horizontal','Position',[0.5 0.5 0 0],'Box','off')
% 
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% %% All DC
% 
% tooHigh = find(initialcond(end,2,:)>6);
% tooLow = find(initialcond(end,2,:)<-2);
% initialcond_noExtremes = initialcond;
% initialcond_noExtremes(:,:,[tooLow;tooHigh]) = [];
% 
% hFig = figure;
% scale = 3;
% numColumns = n;
% numRows = 1;
% set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
% 
% i = 2;
% plot(plotDates,squeeze(initialcond_noExtremes(:,i,:))), hold on;
% plot(plotDates,Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;
% % ylabel(priorTitle)
% % title('Inflation')
% set(gca,'FontSize',18)
% 
% %% DC top 10
% 
% hFig = figure;
% scale = 3;
% numColumns = n;
% numRows = 1;
% set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])
% VARnames = {'GDP growth','Inflation','Fed funds rate'};
% i = 2;
% % subplot(1,n,i)
% plot(plotDates,squeeze(initialcond(:,i,solutions))), hold on;
% plot(plotDates,Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;
% % title(VARnames{i})
% set(gca,'FontSize',18)
% 
% 
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% filename = strcat('Figures/DC_FP10_',priorName);
% print(hFig, filename, '-dpdf','-r0')
% 
% 
% 
