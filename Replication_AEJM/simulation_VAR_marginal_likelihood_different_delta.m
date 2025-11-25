%%% What drives the recent surge in inflation? The historical decomposition roller
%%% Replication codes
%%% Figure 5.(c)

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

%% Generate data
clear;
rng(10)

% Simulate data for a bi-variate VAR(1):

T = 150;

n=2; % Number of variables 

A=[0.95 -0.3;0.3 0.4]; % Autoregressive matrix A
C=[0.4;0.5];
B0=[0.6 0;-0.3 0.5]; % Structural matrix B0

e=standardize(randn(T,n)); % Generate structural errors N(0,1)

x(:,1)=zeros(n,1)+C; % Data

u(:,1)=zeros(n,1); % Reduced-form errors

for t=2:T
    u(:,t)=B0*e(t,:)';
    x(:,t)=C+A*x(:,t-1)+u(:,t); 
end

data = x';
u = u';

%% Estimate VAR model using OLS
[T, n] = size(data);

sur=1;
% Prepare lagged variables
p = 1; % Number of lags

c=1;

delta_all = 0.0001:0.0001:0.4;

y0_bar_factor = [1,1.05];

y0_bar_all = [mean(data,1);...
              mean(data,1)*1.05];

ML_all = nan(size(delta_all,2),size(y0_bar_factor,2));

for yb_i = 1:size(y0_bar_all,1)
    y0_bar = y0_bar_all(yb_i,:);
    for d_i = 1:size(delta_all,2)
        
        if sur==1
            [Y,X,Y_initial] = SUR_sur(data,p,c,y0_bar,delta_all(d_i));
        else
            Y = data;
            lagged_Y = lagmatrix(Y, 1:p); % Create lagged data
            lagged_Y = lagged_Y(p+1:end, :); % Remove NaNs from the beginning
            Y = Y(p+1:end, :); % Align Y
            
        %     Design matrix for regression (lagged variables stacked)
            X = reshape(lagged_Y, [], n * p);
            
            % Add a constant term (intercept)
            X = [ones(size(X, 1), 1), X];                                 
        end
        
        % Estimate VAR coefficients using OLS
        B_hat = (X' * X) \ (X' * Y); % (k x n) coefficient matrix

        % Calculate residuals
        residuals = Y - X * B_hat;
        
        % Compute covariance matrix of residuals
        Sigma_hat = (residuals' * residuals) / (T - n * p - 1);
        
        
        % Marginal Likelihood Calculation
        % Set up prior parameters
        prior_mean = zeros(n * p + 1, n); % Mean of prior for coefficients
        prior_precision = eye(n * p + 1); % Precision matrix of prior for coefficients
        prior_df = n + 1; % Prior degrees of freedom
        prior_S = eye(n); % Prior scale matrix for residual covariance
        
        
        % Posterior parameters
        posterior_precision = prior_precision + X' * X;
        posterior_mean = posterior_precision \ (prior_precision * prior_mean + X' * Y);
        posterior_df = prior_df + T;
        
        posterior_S = prior_S + residuals' * residuals + ...
                      (prior_mean - B_hat)' * prior_precision * (prior_mean - B_hat);
        
        
        % Log marginal likelihood
        log_det_prior_precision = log(det(prior_precision));
        log_det_posterior_precision = log(det(posterior_precision));
        log_det_Sigma_hat = log(det(Sigma_hat));
        log_det_prior_S = log(det(prior_S));
        log_det_posterior_S = log(det(posterior_S));
        
        log_marginal_likelihood = ...
            - (T / 2) * log(2 * pi) ...
            - 0.5 * log_det_prior_precision ...
            + 0.5 * log_det_posterior_precision ...
            - (posterior_df / 2) * log_det_posterior_S ...
            + (prior_df / 2) * log_det_prior_S ...
            - 0.5 * T * log_det_Sigma_hat;
        

        ML_all(d_i,yb_i) = log_marginal_likelihood;
    end
end

%% Figure 5.(c): Marginal likelihood
hFig = figure;
for ii = 1:size(y0_bar_factor,2)
    plot(delta_all,ML_all(:,ii)); hold on
end

lgd = legend({'0%','5%'},'Box','off','FontSize', 18,'Orientation','horizontal','Location','northeast');

title(lgd, ['$\bar{Y}_0$',', deviation from sample mean:'],'Interpreter', 'latex');
    set(gca, 'FontSize', 18, 'Box', 'off', 'Layer', 'top', 'LineWidth', 2, ...
             'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 

ylabel('log marginal likelihood')
xlabel('\delta')

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filename = 'Figures\Fig5c';
print(hFig, filename, '-dpdf','-r0')

