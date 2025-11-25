function Xcond = VARcf_DKcksV2(X,p,beta,Su,nDraws,LinComb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes conditional forecasts for the missing observations in X using a
% VAR and Kalman filter and smoother
%
% INPUTS
% X - matrix of observable variables
% p - # lags in VAR
% beta - coefficients of the VAR
% Su - covariance matrix of the VAR
% nDraws - (optional) number of draws - if == 0 then simple Kalman smoother
%           is run, otherwise nDraws draws of the states are done
%
% OUTPUTS
% Xcond - nans from X are replaced by the conditional forecasts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5
   nDraws = 0;
end

%% Data

N = size(Su,2);

T = size(X,1);

if nargin==6
    q = size(LinComb,2);
    CCadd = zeros(q,N*p); CCadd(:,1:N) = LinComb';
else
    q = 0;
    Cadd = [];
end;


idxNaN = any(isnan(X),2);
idxNaNcs = cumsum(idxNaN(end:-1:1));
nNaNs = sum(idxNaNcs==(1:T)');

% end part with missing observations
Xub = X(end-nNaNs:end,:);
% balanced part
X = X(1:end-nNaNs,:);
Xinit = X(:,q+1:end);

%% State space representation

% transition equation
AA = zeros(N*p);
AA(1:N,1:N*p) = beta(1:end-1,:)'; % autoregressive coefficients
AA(N+1:N*p,1:N*(p-1)) = eye(N*(p-1));
c2 = [beta(end,:)';zeros(N*(p-1),1)];% constant 

% measurement equation
CC = zeros(N,N*p); CC(:,1:N) = eye(N); CC = [CCadd; CC];
QQ = zeros(N*p); QQ(1:N,1:N) = Su;
c1 = zeros(N+q,1);% constant



% initialisation of the Kalman filter
initx = lagmatrix(Xinit,0:p-1);
initx = initx(end,:)';
initV = eye(length(initx))*1e-7;


%% Conditional forecasts

% to save computation time Kalman filter will be run only on the part of data
% with missing observations
yinput = Xub(2:end,:);

Tub = size(yinput,1);

if nDraws == 0 % point forecast
    % Kalman filter and smoother
    [xsmooth] = runKF_DK(yinput', AA, CC, QQ, diag(ones(1,N+q)*1e-12), initx, initV,c1,c2);

    Xcond = [Xinit;xsmooth(1:N,:)']*CC(:,1:N)';
else
    Xcond = nan(T,N,nDraws);
    % Durbin and Koopman simulation smoother
    for kg =1:nDraws
        aplus = nan(N*p,Tub);
        yplus = nan(N,Tub);
        for t = 1:Tub
            aplus(:,t) =  AA*initx+...
                [mvnrnd(zeros(N,1),Su,1)';zeros(N*(p-1),1)]+c2;
            initx = aplus(:,t);
            %         [y_t,Z_t,G_t,M_t] = MissData(yinput(t,:)',CC,G);
            yplus(:,t) = CC*aplus(:,t)+c1;
            %         yplus(isnan(yinput(t,:)),t) = nan;

        end
        ystar = yinput'-yplus;
        ahatstar = runKF_DK(ystar, AA, CC, QQ, diag(ones(1,N+q)*1e-12), zeros(size(initx)), initV,zeros(N,1),zeros(size(initx)));
        atilda = ahatstar+aplus;
        Xcond(:,:,kg) = [Xinit;atilda(1:N,:)]'*CC(:,1:N)';
    end
end