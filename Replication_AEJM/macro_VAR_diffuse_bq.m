%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figures 3.(a) and A-1.(a)

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed

saveFolder = 'Figures';
priorName = 'diffuse_bq';

%% Load data

makedata; % organizes the raw data into a T x n matrix and stores dates

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

hor = 17; % IRF length
candidateirf=zeros(n,n,hor,draws); % candidate impulse response
candidateirf_cum=zeros(n,n,hor,draws); % candidate impulse response

% Set up 4-D matrices for IRFs to be filled in the loop:

C=zeros(n,n,hor,draws); 
D=zeros(n,n,hor,draws);

% Structural shocks:

eta=zeros(T-p,n,draws);

% IRF:

for k=1:draws

for j=1:hor
    BigC=BigA(:,:,k)^(j-1);
    C(:,:,j,k)=BigC(1:n,1:n); % IRFs of the Wold representation
end

% Cholesky factorization:

BigLR=inv(eye(n*p,n*p)-BigA(:,:,k));
C1=BigLR(1:n,1:n);

S=chol(C1*Sigma(:,:,k)*C1','lower'); % lower triangular matrix
K=C1\S;

for i=1:hor
    candidateirf(:,:,i,k)=C(:,:,i,k)*K; % Cholesky IRFs
end

candidateirf_cum(:,:,:,k)=cumsum(candidateirf(:,:,:,k),3);

end

% Compute historical decompositions
hist_decomp;

% Obtain the draws closest to the point-wise median IRF
fry_pagan;


%% Figures
plotDates = dates(p+1:end,:);

%% Figure 3.(a): All deterministic components (DC)

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
%title('All posterior draws')

set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% Save figure

filename = strcat(saveFolder,'/Fig3a');
print(hFig, filename, '-dpdf','-r0')

%% Figure A-1.(a): IRF, Blanchard and Quah

conf=68;
colorBNDS=[0.7 0.7 0.7];
VARnames={'GDP growth'; 'Inflation'};
Shocknames={'Aggregate supply';'Aggregate demand'};


% 3 closest: Max/min
hFig = figure;
scale = 4;
numColumns = n;
numRows = n;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])

for k=1:n % variable k
    for j=1:n % shock j
        
    subplot(n,n,j+n*k-n)
    fill([0:hor-1 fliplr(0:hor-1)]' ,[reshape(permute(max(candidateirf(k,j,:,:),[],4),[3 2 1]),hor,1,[]); flipud(reshape(permute(min(candidateirf(k,j,:,:),[],4),[3 2 1]),hor,1,[]))],...
        [.93 .93 .93],'EdgeColor','None'); hold on;
    fill([0:hor-1 fliplr(0:hor-1)]' ,[reshape(permute(prctile(candidateirf(k,j,:,:),(100+conf)/2,4),[3 2 1]),hor,1,[]); flipud(reshape(permute(prctile(candidateirf(k,j,:,:),(100-conf)/2,4),[3 2 1]),hor,1,[]))],...
        [.6 .6 .6],'EdgeColor','None'); hold on;
%     plot(0:hor-1,reshape(permute(prctile(candidateirf(k,j,:,:),(100-conf)/2,4),[3 2 1]),hor,1,[]),'LineWidth',1.5,'Color','k'); hold on;
    plot(0:hor-1,reshape(permute(prctile(candidateirf(k,j,:,:),50,4),[3 2 1]),hor,1,[]),'LineWidth',3.5,'Color','k'); hold on;
    plot(0:hor-1,candidateirf_wold_fp(:,j+n*k-n),'LineWidth',3.5,'Color','r'); hold on;
    plot(0:hor-1,candidateirf_wold_fp_2(:,j+n*k-n),'LineWidth',3.5,'Color','g'); hold on;
    plot(0:hor-1,candidateirf_wold_fp_3(:,j+n*k-n),'LineWidth',3.5,'Color','b'); hold on;
%     plot(0:hor-1,reshape(permute(prctile(candidateirf(k,j,:,:),(100+conf)/2,4),[3 2 1]),hor,1,[]),'LineWidth',1.5,'Color','k'); hold on;
    line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5); hold off;
        set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
    
    if j ==1
        ylabel(strcat(VARnames{k}), 'FontSize', 18);
    end
    
    if k == 1
        title(Shocknames{j}, 'FontSize', 18);
    end
 
    xlim([0 hor-1]);
    
    end
end
legend({'min/max','68% credible set','Pointwise median','Fry-Pagan draw','2nd closest','3rd closest'},'FontSize',16,'Orientation','Horizontal','Position',[0.5 0.5 0 0],'Box','off')

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/FigA1a');
print(hFig, filename, '-dpdf','-r0')



