%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figure 5.(d)

clear; close all; clc;

MTBfolder = 'Functions/';

addpath(strcat(MTBfolder,'getfreddata-matlab-master'))
addpath(strcat(MTBfolder,'Functions_Primiceri/GLPreplicationWeb'))
addpath(strcat(MTBfolder,'Functions_Primiceri/BGLreplication'))
addpath(strcat(MTBfolder,'Functions_Primiceri/GLPreplicationWeb/subroutines'))
addpath(strcat(MTBfolder,'canova/BVAR_-master/bvartools'))
addpath(strcat(MTBfolder,'matlab_functions'))
addpath(strcat(MTBfolder,'binning'))

set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

%% Simulate data for a VAR(1) with three variables:
clear;
rng(10)

T=250; % Length of the data (the 100 first observations are used as pre-sample)


n=2; % variables 

A=[0.9 -0.3;0.3 0.4]; % Autoregressive matrix
C=[0.4;0.5];
B0=[0.6 0;-0.3 0.5]; % Structural matrix B0

e=standardize(randn(T,n)); % Generate structural errors N(0,1)

x(:,1)=zeros(n,1)+C; % Data
%x(:,1)=zeros(n,1); % Data
u(:,1)=zeros(n,1); % Reduced-form errors

for t=2:T
        
        u(:,t)=B0*e(t,:)';
        x(:,t)=C+A*x(:,t-1)+u(:,t);
        
end

x = x';
u = u';

uncond_mean = (eye(n)-A)^(-1)*C;
% mean(x,1)
pre_vardata = x(1:100,:);
x = x(101:end,:);

%% Estimate a VAR(1) with three variables using simulated data:

code = 1;

vardata = x;

minnesota = 1;
sur = 0;
noc = 0;

lambda = []; %0.2;
theta  = []; % 0.0001;
mu     = [];%1;

iid = [1,2];

p=1;
c=1;



delta_all = 0.001:0.005:0.2;

disp_table = nan(size(delta_all,2),2,5);

for y0_i = 1:5
    if y0_i==1
        y0_bar = mean(vardata,1);
    elseif y0_i==2
        y0_bar = uncond_mean';
    elseif y0_i==3
        y0_bar = mean(vardata(1:p,:),1);
    elseif y0_i==4
        y0_bar = mean(pre_vardata,1);
    elseif y0_i==5
        y0_bar = [0,0];
    end
    T = size(vardata,1);
    % Set up the loop for each draw :
    T = T+1;
    draws=1000;

    PI=zeros(n*p+c,n,draws);
    BigA=zeros(n*p,n*p,draws);
    Sigma=zeros(n,n,draws);
    errornorm=zeros(T-p,n,draws);
    fittednorm=zeros(T-p,n,draws);
    delta_counter = 0;
    for d_i = 1:size(delta_all,2)
        delta = delta_all(d_i);
        for i=1:draws    
            stable=-1;    
            while stable<0
    
            [PI(:,:,i),BigA(:,:,i),Sigma(:,:,i),errornorm(:,:,i),fittednorm(:,:,i)]=BVAR_sur(vardata,p,c,y0_bar,delta);
    
            if abs(eig(BigA(:,:,i)))<1
                stable=1; % keep only stable draws
            end  
            end   
        end
   
    
        %% Historical decomposition and initial conditions:
        
        % Data:
        
        Y=x(p+1:end,:); % Discard first p observations
        
        % Stochastic component:
        T = T-1;
        err=errornorm; % Reduced-form residuals
        HH=zeros(T-p,n,n,draws); % Store the impact IRFs ^t 
        At=zeros(T-p,n,n,draws); % Store A^t 
        hist=zeros(T-p,n,n,draws);
        
        for k=1:draws
        for t=1:T-p  
            
                BigC=(BigA(:,:,k))^(t-1);
                BigH=BigC(1:n,1:n);
                At(t,:,:,k)=BigH;
                
                for j=1:n
                    
                HH(t,:,j,k)=BigH(j,1:n);
                
                end       
        end
        
        %HH=cumsum(HH,1); % If variables are in first differences
        
        for t=1:T-p
            for i=1:n
                for j=1:n
                
                  hist(t,i,j,k)=HH(1:t,i,j,k)'*flipud(err(1:t,i,k));
                  
                end        
            end   
        end
        end
        
        % Deterministic component:
        
        initialcond=zeros(T-p,n,draws);
        Aty0=zeros(T-p,n,draws);
        
        for k=1:draws % Draws
        for t=1:T-p % Time
        
            for i=1:n % Variables 
                
                initialcond(t,i,k)=Y(t,i)-sum(hist(t,:,i,k),2);  
                
            end   
        
            Aty0(t,:,k)=squeeze(At(t,:,:,k))*x(1,:)';
        
        end
        end
        
        % save('simul_150_02_DC_all_min','initialcond','Y')
        
        %%
        dc_var1 = squeeze(initialcond(:,1,:));
        std_dc_var1 = std(dc_var1(end,:));
        std_var1 = std(vardata(:,1));
        
        dispersion = std_dc_var1/std_var1;
        
        disp_table(d_i,1,y0_i) = delta;
        disp_table(d_i,2,y0_i) = dispersion;
    end
end




%% Fig 5.(d): Dispersion of deterministic components
hFig = figure;
plot(disp_table(:,1,1),disp_table(:,2,1)); hold on
plot(disp_table(:,1,2),disp_table(:,2,2))
plot(disp_table(:,1,3),disp_table(:,2,3))
plot(disp_table(:,1,4),disp_table(:,2,4))
plot(disp_table(:,1,5),disp_table(:,2,5))

legend({'Sample mean (0.56)','Steady state mean (0.60) ','Average of first $p$ observations (-2.24)','Zero','Presample mean (0.68)'},'Box','off','Location','southeast','FontSize',18,'Interpreter','latex')
% ylabel('Relative dispersion of DC $\frac{\sigma^{DC}_{1,T}}{\sigma_{1}}$','Interpreter','latex')
ylabel('Relative dispersion of $DC_T$','Interpreter','latex')
xlabel('\delta')
    set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = 'Figures\Fig5d';
print(hFig, filename, '-dpdf','-r0')
