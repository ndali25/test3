%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figures 3.(f) and E-1

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed

MTBfolder = 'C:\Users\06132027\Nord universitet Dropbox\Pål Ulvedal\MatlabTB\';

addpath(strcat(MTBfolder,'getfreddata-matlab-master'))
addpath(strcat(MTBfolder,'Functions_Primiceri\GLPreplicationWeb'))
addpath(strcat(MTBfolder,'Functions_Primiceri\BGLreplication'))
addpath(strcat(MTBfolder,'Functions_Primiceri\GLPreplicationWeb\subroutines'))
addpath(strcat(MTBfolder,'canova\BVAR_-master\bvartools'))
addpath(strcat(MTBfolder,'arias_sign_zero'))
addpath(strcat(MTBfolder,'matlab_functions'))
addpath(strcat(MTBfolder,'binning'))

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
draws = 10000; % Number of final draws
ndDisc = 0; % Number of draws to discard


saveFolder = 'Figures';
priorName = 'diffuse_large_model';

%% SVAR with sign restrictions - diffuse prior with large system

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
% Set restricions

noConstant = 0;
% Idenfify shocks, reduced-form

VARnames={'GDP'; 'Inflation';'Investment';'Fed funds rate';'Real wage'};
Shocknames={'Other demand';'Investment';'Monetary policy';'Labor supply';'Other supply'};

% Compute historical decompositions

Y=data(p+1:end,:); 
err=errornorm; % reduced form residuals
xi=zeros(size(errornorm,1),size(errornorm,2),size(errornorm,3));
HH=zeros(T-p,n,n,draws); % store the impact IRFs ^t 
histdec=zeros(T-p,n,n,draws);

for k=1:draws % draw k
    
    xi(:,:,k)=err(:,:,k); 
    
    for t=1:T-p
            
            BigC=(BigA(:,:,k))^(t-1);
            BigH=BigC(1:n,1:n); % C * S * H
            
            for j=1:n % variable j
                
            HH(t,:,j,k)=BigH(j,1:n);
            
            end
            
    end
    
    for t=1:T-p
       
        for i=1:n % shock i 
            
            for j=1:n % variable j
            
              histdec(t,i,j,k)=HH(1:t,i,j,k)'*flipud(xi(1:t,i,k));
              
            end
            
        end
        
    end

end

% Deterministic component:

initialcond=zeros(T-p,n,draws);

for k=1:draws % draws
    for t=1:T-p % time
        for i=1:n % variables 
            
            initialcond(t,i,k)=Y(t,i)-sum(histdec(t,:,i,k),2);  
            
        end   
    end
end


%% Figures
plotDates = datetime(dates(p+1:end,:));


%% Fig 3.(f): All DC for inflation

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

filename = strcat(saveFolder,'/Fig3f');
print(hFig, filename, '-dpdf','-r0')


%% Fig E-1: All DCs all variables

tooHigh = find(initialcond(end,2,:)>2);
tooLow = find(initialcond(end,2,:)<-0.5);
initialcond_noExtremes = initialcond;
initialcond_noExtremes(:,:,[tooLow;tooHigh]) = [];

hFig = figure;
scale = 5;
numColumns = 2;
numRows = 3;
set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])


for i = 1:5
    subplot(numRows,numColumns,i)
    plot(plotDates,squeeze(initialcond_noExtremes(:,i,:))), hold on;
    plot(plotDates,Y(:,i),'-','LineWidth',2,'Color','k'), axis tight;
%     ylim([-0.5, 2.5])
    set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
    title(VARnames{i})
end

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = strcat(saveFolder,'/FigE1');
print(hFig, filename, '-dpdf','-r0')

