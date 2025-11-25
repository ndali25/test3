%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figures 3.(c) and 10.(b)

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed

saveFolder = 'Figures';
priorName = 'NIW';

%% Load data

makedata; % organizes the raw data into a T x n matrix and stores dates

%% SVAR with sign restrictions - NIW prior

draws = 10000; % Number of final draws
ndDisc = 3000; % Number of draws to discard

p = 4; % Lags

options.priors.name = 'Conjugate';
options.K = draws;
BVAR = bvar_(data,p,options);
Bdraw_all        = [BVAR.Phi_draws(end,:,:); BVAR.Phi_draws(1:end-1,:,:)];
Sigmadraw_all     = BVAR.Sigma_draws;


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

%% Figure 3.(c): All DC

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

% 
filename = strcat(saveFolder,'/Fig3c');
print(hFig, filename, '-dpdf','-r0')


%% Figure 10.(b): Median historical decomposition - 2019 onwards

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
filename = strcat(saveFolder,'/Fig10b');
print(hFig, filename, '-dpdf','-r0')

%% Uncomment below for IRFs and HDs for NIW:

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
%     set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
%              'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
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
% set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
%              'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
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
% set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
%              'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
% 
% 
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% filename = strcat('Figures/DC_FP10_',priorName);
% print(hFig, filename, '-dpdf','-r0')
% 
% %% HD top 3
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
% % if k==1
% %     ylabel(priorTitle)
% % end
% title(titles{k})
% 
% % ylim(yLimit)
% pos = get(gca,'Position');
% pos(2) = pos(2)+.1;
% pos(4) = pos(4)-.2;
% set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
%              'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
% end
% legend({'Deterministic component','Supply','Demand'},'Orientation','Horizontal','Position',[0.5 0.02 0 0],'Box','off')
% 
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% 
