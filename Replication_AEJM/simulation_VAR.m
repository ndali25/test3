%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figures 4, C-2 and C-3

clc; clear;

%% General settings

addpath(genpath('Functions'))
set(0, 'defaultAxesFontName', 'Times'); % Clean and modern font
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5);

rng(1) % Set seed

saveFolder = 'Figures';
priorName = 'diffuse_cholesky';


%% Simulate data for a VAR(1) with three variables:
rng(12345)

T_large=500; % Length of the data
T_med=150; % Length of the data
T_small=80; % Length of the data
% T = 180;


n=2; % variables 

%A=[0.6 -0.3;0.3 0.4]; % Autoregressive matrix A
A=[0.95 -0.3;0.3 0.4]; % Autoregressive matrix A
A2=[0.6 -0.3;0.3 0.4]; % Autoregressive matrix A
C=[0.4;0.5];
B0=[0.6 0;-0.3 0.5]; % Structural matrix B0

e=standardize(randn(T_large,n)); % Generate structural errors N(0,1)

x=zeros(n,T_large); % Data
x2=zeros(n,T_large); % Data
x(:,1)=C; % Data
x2(:,1)=C; % Data
%x(:,1)=zeros(n,1); % Data
u=zeros(n,T_large); % Reduced-form errors

for t=2:T_large
        
        u(:,t)=B0*e(t,:)';
        x(:,t)=C+A*x(:,t-1)+u(:,t);
        x2(:,t)=C+A2*x2(:,t-1)+u(:,t);
        
end

x = x';
x2 = x2';
u = u';

%% Generate true IRF

hor = 20;

IRF_true = nan(hor+1,n,n);

for ii = 0:hor
     temp_IRF = A^ii*B0;
     for jj = 1:n
         IRF_true(ii+1,:,jj) = temp_IRF(jj,:); 
         IRF_true(ii+1,:,jj) = temp_IRF(jj,:);
     end
end


%% Estimate a VAR(1) (persistent case) with Diffuse prior
p=1;
c=1;
noConstant = 0;
[f,sr,var_pos] = var_restrictions('cholesky_sim');

% Draws set up:
draws=1000;

% Small sample
estimation_name = 'diffuse_small';
vardata = x(1:T_small,:);
data.(estimation_name) = vardata;

options.priors.name = 'Jeffrey';
options.K = draws;
BVAR = bvar_(vardata,p,options);
Bdraw_all        = [BVAR.Phi_draws(end,:,:); BVAR.Phi_draws(1:end-1,:,:)];
Sigmadraw_all     = BVAR.Sigma_draws;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

% Medium sample
estimation_name = 'diffuse_med';
vardata = x(1:T_med,:);
data.(estimation_name) = vardata;

options.priors.name = 'Jeffrey';
options.K = draws;
BVAR = bvar_(vardata,p,options);
Bdraw_all        = [BVAR.Phi_draws(end,:,:); BVAR.Phi_draws(1:end-1,:,:)];
Sigmadraw_all     = BVAR.Sigma_draws;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

% Large sample
estimation_name = 'diffuse_large';
vardata = x(1:T_large,:);
data.(estimation_name) = vardata;

options.priors.name = 'Jeffrey';
options.K = draws;
BVAR = bvar_(vardata,p,options);
Bdraw_all        = [BVAR.Phi_draws(end,:,:); BVAR.Phi_draws(1:end-1,:,:)];
Sigmadraw_all     = BVAR.Sigma_draws;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));


%% Estimate a VAR(1) (less persistent case) with Diffuse prior
p=1;
c=1;
noConstant = 0;
[f,sr,var_pos] = var_restrictions('cholesky_sim');

% Draws set up:
draws=1000;

% Small sample
estimation_name = 'diffuse_small_lp';
vardata = x2(1:T_small,:);
data.(estimation_name) = vardata;

options.priors.name = 'Jeffrey';
options.K = draws;
BVAR = bvar_(vardata,p,options);
Bdraw_all        = [BVAR.Phi_draws(end,:,:); BVAR.Phi_draws(1:end-1,:,:)];
Sigmadraw_all     = BVAR.Sigma_draws;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

% Medium sample
estimation_name = 'diffuse_med_lp';
vardata = x2(1:T_med,:);
data.(estimation_name) = vardata;

options.priors.name = 'Jeffrey';
options.K = draws;
BVAR = bvar_(vardata,p,options);
Bdraw_all        = [BVAR.Phi_draws(end,:,:); BVAR.Phi_draws(1:end-1,:,:)];
Sigmadraw_all     = BVAR.Sigma_draws;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

% Large sample
estimation_name = 'diffuse_large_lp';
vardata = x2(1:T_large,:);
data.(estimation_name) = vardata;

options.priors.name = 'Jeffrey';
options.K = draws;
BVAR = bvar_(vardata,p,options);
Bdraw_all        = [BVAR.Phi_draws(end,:,:); BVAR.Phi_draws(1:end-1,:,:)];
Sigmadraw_all     = BVAR.Sigma_draws;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));



%% Estimate a VAR(1) (persistent case) with Minnesota
sur = 0;
noc = 0;
noConstant = 0;

iid = [1,2];

p=1;
c=1;
[f,sr,var_pos] = var_restrictions('cholesky_sim');

% Draws set up:
draws=1000;
ndDisc = 1000;

% Small sample
estimation_name = 'minnesota_small';
vardata = x(1:T_small,:);
data.(estimation_name) = vardata;
y0_bar = mean(vardata,1);

res.(estimation_name)    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
Bdraw_all = res.(estimation_name).mcmc.beta;
Sigmadraw_all = res.(estimation_name).mcmc.sigma;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

% Medium sample
estimation_name = 'minnesota_med';
vardata = x(1:T_med,:);
data.(estimation_name) = vardata;
y0_bar = mean(vardata,1);

res.(estimation_name)    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
Bdraw_all = res.(estimation_name).mcmc.beta;
Sigmadraw_all = res.(estimation_name).mcmc.sigma;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

% Large sample
estimation_name = 'minnesota_large';
vardata = x(1:T_large,:);
data.(estimation_name) = vardata;
y0_bar = mean(vardata,1);

res.(estimation_name)    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
Bdraw_all = res.(estimation_name).mcmc.beta;
Sigmadraw_all = res.(estimation_name).mcmc.sigma;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

%% Estimate a VAR(1) (persistent case) with Minnesota and SUR
sur = 1;
noc = 0;
noConstant = 0;

iid = [1,2];

p=1;
c=1;
[f,sr,var_pos] = var_restrictions('cholesky_sim');

% Draws set up:
draws=1000;
ndDisc = 1000;

% Small sample
estimation_name = 'sur_small';
vardata = x(1:T_small,:);
data.(estimation_name) = vardata;
y0_bar = mean(vardata,1);

res.(estimation_name)    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
Bdraw_all = res.(estimation_name).mcmc.beta;
Sigmadraw_all = res.(estimation_name).mcmc.sigma;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

% Medium sample
estimation_name = 'sur_med';
vardata = x(1:T_med,:);
data.(estimation_name) = vardata;
y0_bar = mean(vardata,1);

res.(estimation_name)    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
Bdraw_all = res.(estimation_name).mcmc.beta;
Sigmadraw_all = res.(estimation_name).mcmc.sigma;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

% Large sample
estimation_name = 'sur_large';
vardata = x(1:T_large,:);
data.(estimation_name) = vardata;
y0_bar = mean(vardata,1);

res.(estimation_name)    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
Bdraw_all = res.(estimation_name).mcmc.beta;
Sigmadraw_all = res.(estimation_name).mcmc.sigma;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

%% Estimate a VAR(1) (less persistent case) with Minnesota and SUR
sur = 1;
noc = 0;
noConstant = 0;

iid = [1,2];

p=1;
c=1;
[f,sr,var_pos] = var_restrictions('cholesky_sim');

% Draws set up:
draws=1000;
ndDisc = 1000;

% Small sample
estimation_name = 'sur_small_lp';
vardata = x2(1:T_small,:);
data.(estimation_name) = vardata;
y0_bar = mean(vardata,1);

res.(estimation_name)    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
Bdraw_all = res.(estimation_name).mcmc.beta;
Sigmadraw_all = res.(estimation_name).mcmc.sigma;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

% Medium sample
estimation_name = 'sur_med_lp';
vardata = x2(1:T_med,:);
data.(estimation_name) = vardata;
y0_bar = mean(vardata,1);

res.(estimation_name)    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
Bdraw_all = res.(estimation_name).mcmc.beta;
Sigmadraw_all = res.(estimation_name).mcmc.sigma;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));

% Large sample
estimation_name = 'sur_large_lp';
vardata = x2(1:T_large,:);
data.(estimation_name) = vardata;
y0_bar = mean(vardata,1);

res.(estimation_name)    = bvarGLP_y0(vardata, p, 'noc',noc,'sur', sur, 'mcmc', 1, 'MCMCconst', 2, 'MNalpha',0,'MNpsi',0,'pos', iid, 'Ndraws', draws+ndDisc, 'Ndrawsdiscard', ndDisc, 'hyperpriors', 0,'y0_custom',y0_bar); % ,'pos', [1,2,3,4]
Bdraw_all = res.(estimation_name).mcmc.beta;
Sigmadraw_all = res.(estimation_name).mcmc.sigma;

R.(estimation_name) = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos);
[histdec.(estimation_name), det_comp.(estimation_name)] = get_hist_decomp(vardata,p,c,Bdraw_all,R.(estimation_name));


%% Figure 4.(a): Estimated deterministic components of Y1,t. Diffuse prior

% variables = [2,1];
estimations = {'diffuse_small','diffuse_med','diffuse_large'};
estimations_lp = {'diffuse_small_lp','diffuse_med_lp','diffuse_large_lp'};
titles = {'T = 80', 'T = 150', 'T = 500'};
% ylabels = {'More persistent','Less persistent'};


hFig = figure;
scale = 5;
numColumns = 3;
numRows = 2;

set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])

kk = 0;
i = 1;
for j = 1:length(estimations_lp)
    kk = kk+1;
    subplot(2,3,kk)
    plot(squeeze(det_comp.(estimations_lp{j})(:,i,:))), hold on;
    plot(data.(estimations_lp{j})(p+1:end,i),'-','LineWidth',2,'Color','k'), axis tight;
    title(titles{j})
    if j==1
        ylabel('Less persistent')
    end
    ylim([-3, 5])
    set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
end
for j = 1:length(estimations)
    kk = kk+1;
    subplot(2,3,kk)
    plot(squeeze(det_comp.(estimations{j})(:,i,:))), hold on;
    plot(data.(estimations{j})(p+1:end,i),'-','LineWidth',2,'Color','k'), axis tight;
    title(titles{j})
    if j==1
        ylabel('More persistent')
    end
    ylim([-3, 5])
    set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
end


set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = 'Figures/Fig4a';
print(hFig, filename, '-dpdf','-r0')

%% Figure 4.(b): Estimated deterministic components of Y1,t. SUR prior

% variables = [2,1];
estimations = {'sur_small','sur_med','sur_large'};
estimations_lp = {'sur_small_lp','sur_med_lp','sur_large_lp'};
titles = {'T = 80', 'T = 150', 'T = 500'};
% ylabels = {'More persistent','Less persistent'};


hFig = figure;
scale = 5;
numColumns = 3;
numRows = 2;

set(hFig,'position',[200 200 160*scale*numColumns 90*scale*numRows])

kk = 0;
i = 1;
for j = 1:length(estimations_lp)
    kk = kk+1;
    subplot(2,3,kk)
    plot(squeeze(det_comp.(estimations_lp{j})(:,i,:))), hold on;
    plot(data.(estimations_lp{j})(p+1:end,i),'-','LineWidth',2,'Color','k'), axis tight;
    title(titles{j})
    if j==1
        ylabel('Less persistent')
    end
    ylim([-3, 5])
    set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
end
for j = 1:length(estimations)
    kk = kk+1;
    subplot(2,3,kk)
    plot(squeeze(det_comp.(estimations{j})(:,i,:))), hold on;
    plot(data.(estimations{j})(p+1:end,i),'-','LineWidth',2,'Color','k'), axis tight;
    title(titles{j})
    if j==1
        ylabel('More persistent')
    end
    ylim([-3, 5])
    set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
end


set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


filename = 'Figures/Fig4b';
print(hFig, filename, '-dpdf','-r0')
%% Fig C-2: Posterior distributions of impulse responses, Minnesota and Minnesota and SUR priors, simulated data

R_min = R.minnesota_med;
R_min_sur = R.sur_med;

x_lab = 0:20;

band = 0.95;
% Minnesota IRF
R_min_short = R_min(:,1:hor+1,:,:);
IRF_min = squeeze(quantile(R_min_short,0.5,3));
IRF_min_low = squeeze(quantile(R_min_short,0.5-band/2,3));
IRF_min_high = squeeze(quantile(R_min_short,0.5+band/2,3));

% Minnesota-SUR IRF
R_min_sur_short = R_min_sur(:,1:hor+1,:,:);
IRF_min_sur = squeeze(quantile(R_min_sur_short,0.5,3));
IRF_min_sur_low = squeeze(quantile(R_min_sur_short,0.5-band/2,3));
IRF_min_sur_high = squeeze(quantile(R_min_sur_short,0.5+band/2,3));

hFig = figure;
subplot(2,2,1)
plot(x_lab,IRF_true(:,1,1),'k'); hold on
plot(x_lab,IRF_min(1,:,1),'r');
plot(x_lab,IRF_min_high(1,:,1),'LineStyle','--','Color','r');
plot(x_lab,IRF_min_low(1,:,1),'LineStyle','--','Color','r');

plot(x_lab,IRF_min_sur(1,:,1),'b');
plot(x_lab,IRF_min_sur_high(1,:,1),'LineStyle','--','Color','b');
plot(x_lab,IRF_min_sur_low(1,:,1),'LineStyle','--','Color','b');
title('Shock 1')
ylabel('Variable 1')
yline(0)
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

subplot(2,2,2)

plot(x_lab,IRF_true(:,2,1),'k'); hold on
plot(x_lab,IRF_min(1,:,2),'r');
plot(x_lab,IRF_min_high(1,:,2),'LineStyle','--','Color','r');
plot(x_lab,IRF_min_low(1,:,2),'LineStyle','--','Color','r');

plot(x_lab,IRF_min_sur(1,:,2),'b');
plot(x_lab,IRF_min_sur_high(1,:,2),'LineStyle','--','Color','b');
plot(x_lab,IRF_min_sur_low(1,:,2),'LineStyle','--','Color','b');
title('Shock 2')
yline(0)
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

subplot(2,2,3)
plot(x_lab,IRF_true(:,1,2),'k'); hold on
plot(x_lab,IRF_min(2,:,1),'r');
plot(x_lab,IRF_min_high(2,:,1),'LineStyle','--','Color','r');
plot(x_lab,IRF_min_low(2,:,1),'LineStyle','--','Color','r');

plot(x_lab,IRF_min_sur(2,:,1),'b');
plot(x_lab,IRF_min_sur_high(2,:,1),'LineStyle','--','Color','b');
plot(x_lab,IRF_min_sur_low(2,:,1),'LineStyle','--','Color','b');
ylabel('Variable 2')
yline(0)
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

subplot(2,2,4)
plot(x_lab,IRF_true(:,2,2),'k'); hold on
plot(x_lab,IRF_min(2,:,2),'r');
plot(x_lab,IRF_min_high(2,:,2),'LineStyle','--','Color','r');
plot(x_lab,IRF_min_low(2,:,2),'LineStyle','--','Color','r');

plot(x_lab,IRF_min_sur(2,:,2),'b');
plot(x_lab,IRF_min_sur_high(2,:,2),'LineStyle','--','Color','b');
plot(x_lab,IRF_min_sur_low(2,:,2),'LineStyle','--','Color','b');
yline(0)
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
legend({'True','Minnesota','95 % band','','Minnesota and SUR','95 % band'},'Box','off','Orientation','horizontal','Position',[0.5 0.02 0 0])


set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print('Figures\FigC2','-dpng','-r0')

%% Fig C-3: Posterior distributions, Flat and SUR priors, simulated data
c_min = res.minnesota_med.mcmc.beta(1,:, :);       % 1x2 constant vector
c_sur = res.sur_med.mcmc.beta(1,:, :);       % 1x2 constant vector
Phi_min = res.minnesota_med.mcmc.beta(2:end,:, :); % 2x2 coefficient matrix
Phi_sur = res.sur_med.mcmc.beta(2:end,:, :); % 2x2 coefficient matrix

I = eye(n);   % 2x2 identity matrix


Y_bar_min = nan(n,draws);
Y_bar_sur = nan(n,draws);
for ii = 1:draws
    Y_bar_min(:,ii) = (I - Phi_min(:,:,ii)') \ c_min(:,:,ii)';
    Y_bar_sur(:,ii) = (I - Phi_sur(:,:,ii)') \ c_sur(:,:,ii)';
end


hFig = figure(Position=[50 50 1500 1500]);

% A-coefficients
subplot(4, 2, 1);
histogram(Phi_min(1,1,:), BinWidth=0.01, FaceColor="blue", FaceAlpha=0.15); hold on;% 30 bins
histogram(Phi_sur(1,1,:), BinWidth=0.01, FaceColor="green", FaceAlpha=0.3); % 30 bins
xlabel('$A_{1,1}$','Interpreter','latex');
ylabel('Frequency');
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

subplot(4, 2, 2);
histogram(Phi_min(1,2,:), BinWidth=0.01, FaceColor="blue", FaceAlpha=0.15); hold on;% 30 bins
histogram(Phi_sur(1,2,:), BinWidth=0.01, FaceColor="green", FaceAlpha=0.3); % 30 bins
xlabel('$A_{1,2}$','Interpreter','latex');
ylabel('Frequency');
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

subplot(4, 2, 3);
histogram(Phi_min(2,1,:), BinWidth=0.01, FaceColor="blue", FaceAlpha=0.15); hold on;% 30 bins
histogram(Phi_sur(2,1,:), BinWidth=0.01, FaceColor="green", FaceAlpha=0.3); % 30 bins
xlabel('$A_{2,1}$','Interpreter','latex');
ylabel('Frequency');
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

subplot(4, 2, 4);
histogram(Phi_min(2,2,:), BinWidth=0.01, FaceColor="blue", FaceAlpha=0.15); hold on;% 30 bins
histogram(Phi_sur(2,2,:), BinWidth=0.01, FaceColor="green", FaceAlpha=0.3); % 30 bins
xlabel('$A_{2,2}$','Interpreter','latex');
ylabel('Frequency');
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

% C-coefficients
subplot(4, 2, 5);
histogram(c_min(1,1,:), BinWidth=0.01, FaceColor="blue", FaceAlpha=0.15); hold on;% 30 bins
histogram(c_sur(1,1,:), BinWidth=0.01, FaceColor="green", FaceAlpha=0.3); % 30 bins
xlabel('$C_{1}$','Interpreter','latex');
ylabel('Frequency');
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

subplot(4, 2, 6);
histogram(c_min(1,2,:), BinWidth=0.01, FaceColor="blue", FaceAlpha=0.15); hold on;% 30 bins
histogram(c_sur(1,2,:), BinWidth=0.01, FaceColor="green", FaceAlpha=0.3); % 30 bins
xlabel('$C_{2}$','Interpreter','latex');
ylabel('Frequency');
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

% (I-A)^-1 C
subplot(4, 2, 7);
histogram(Y_bar_min(1,:), BinWidth=0.01, FaceColor="blue", FaceAlpha=0.15); hold on;% 30 bins
histogram(Y_bar_sur(1,:), BinWidth=0.01, FaceColor="green", FaceAlpha=0.3); % 30 bins
xlabel('$(I-A)^{-1} C$','Interpreter','latex');
ylabel('Frequency');
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

subplot(4, 2, 8);
histogram(Y_bar_min(2,:), BinWidth=0.01, FaceColor="blue", FaceAlpha=0.15); hold on;% 30 bins
histogram(Y_bar_sur(2,:), BinWidth=0.01, FaceColor="green", FaceAlpha=0.3); % 30 bins
xlabel('$(I-A)^{-1} C$','Interpreter','latex');
ylabel('Frequency');
set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 


legend({'Minnesota','Minnesota and Single-unit-root'}, Box="off", Position=[0.5 0.02 0 0], Orientation="horizontal")
% 'Single-unit-root: \delta = 0.1'

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print('Figures\FigC3','-dpdf','-r0')



