function Xcond = DFMcf(X,r,q,p,max_iter,Xest)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes conditional forecasts for the missing observations in X using a
% dynamic factor model and Kalman filter and smoother
% Estimation is based on principal components or ML by EM algorithm
%
% INPUTS
% X - matrix of observable variables
% r - # of static factors
% q - # of dynamic factors
% p - # lags in factor VAR
% max_iter - (optional) #number of EM iterations, if=0 then estimation
%             is based on principal components
% Xest - (optional) matrix for the estimation, default: X
%
% OUTPUTS
% Xcond - nans from X are replace by the conditional forecasts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


OPTS.disp = 0;
[T,N] = size(X);

%% Data
% removing nans at the end of the sample
idxNaN = any(isnan(X),2);
idxNaNcs = cumsum(idxNaN(end:-1:1));
nNaNs = sum(idxNaNcs==(1:T)');

% Unbalanced data
Xub = X;
% Balanced data (for the estimation)
X = X(1:end-nNaNs,:);

% In case different data should be used for the estimation
if nargin>5
    X = Xest;
end

% if max_iter not specified to the PC
if nargin<5
    max_iter = 0;
end

Te = size(X,1);

% standardise data
Mx = mean(X);
Wx = (std(X));
% x is the matrix of standardized observable variables
x = (X-kron(ones(Te,1),Mx))*diag(1./Wx);
xub = (Xub-kron(ones(T,1),Mx))./repmat(Wx,T,1);


%% Estimation
nlag = p-1;

A_temp = zeros(r,r*(nlag + 1))';
I = eye(r*(nlag+1),r*(nlag+1));
A = [A_temp';I(1:end-r,1:end)];

Q = zeros((nlag+1)*r,(nlag+1)*r);
Q(1:r,1:r) = eye(r);

%extract the first r eigenvectors and eigenvalues from cov(x)
[ v, d ] = eigs(cov(x),r,'lm',OPTS);

F = x*v;

beta = (F'*F)\(F'*x);
C = [zeros(N,r) zeros(N,r*(nlag))];
C(:,1:r) = beta';

R = diag(diag(cov(x- F*C(:,1:r)')));


if p > 0    
    z = F;
    Z = [];
    for kk = 1:p
        Z = [Z z(p-kk+1:end-kk,:)]; % stacked regressors (lagged SPC)
    end;
    z = z(p+1:end,:);
    % run the var chi(t) = A*chi(t-1) + e(t);
    A_temp = (Z'*Z)\(Z'*z);        % OLS estimator of the VAR transition matrix
    A(1:r,1:r*p) = A_temp';
    e = z  - Z*A_temp;              % VAR residuals
    H = cov(e);                     % VAR covariance matrix
    if r > 1
        dsp.disp  = 0;
        [D, M]   = eigs(H,q,'lm',dsp);
        D       = D*diag(sign(D(1,:)));
    else
        D = 1;
        M = H;
    end
    Q          = zeros(p*r,p*r);
    Q(1:r,1:r) = D*M*D';
    
end;

initV = reshape(pinv(kron(eye(r*p)-A',eye(r*p)-A'))*Q(:),r*p,r*p);     
initx = zeros(r*p,1);

%% The EM loop
% some auxiliary variables for the iterations
previous_loglik = -inf;
num_iter = 0;
LL = -inf;
converged = 0;
thresh = 1e-4;
decrease = nan(max_iter,1);
LL = nan(max_iter,1);
while (num_iter < max_iter) && ~converged
    [C_new, R_new, A_new, Q_new, initx, initV, loglik] = EMstep(x', A, C, Q, R, initx, initV, r);
    
    C(:,1:r) = C_new;
    R = R_new;
    A = A_new;
    Q = Q_new;

    % Checking convergence
    [converged,decrease(num_iter+1)] = em_converged(loglik, previous_loglik, thresh,1);
     
    LL(num_iter+1) = loglik;
    previous_loglik = loglik;
    num_iter =  num_iter + 1;
end



%% Conditional forecasts
% Kalman filter and smoother
xsmooth = runKF(xub', A, C, Q, R, initx, initV);
% Common component
chi = xsmooth(:,2:end)'*C'*diag(Wx) + kron(ones(T,1),Mx);

% replacing the common component by observations, if available
Xcond = chi;
Xcond(isfinite(Xub)) = Xub(isfinite(Xub)) ;


%--------------------------------------------------------------------------
function  [C_new, R_new, A_new, Q_new, Z_0, V_0, loglik] = EMstep(y, A, C, Q, R, Z_0, V_0, r)

[n,T] = size(y);

% Compute the (expected) sufficient statistics for a single Kalman filter sequence.

%Running the Kalman filter with the current estimates of the parameters
[Zsmooth, Vsmooth, VVsmooth, loglik] = runKF(y, A, C, Q, R, Z_0, V_0);


EZZ = Zsmooth(:,2:end)*Zsmooth(:,2:end)'+sum(Vsmooth(:,:,2:end),3);                        %E(Z'Z)
EZZ_BB = Zsmooth(:,1:end-1)*Zsmooth(:,1:end-1)'+sum(Vsmooth(:,:,1:end-1),3); %E(Z(-1)'Z_(-1))
EZZ_FB = Zsmooth(:,2:end)*Zsmooth(:,1:end-1)'+sum(VVsmooth,3);%E(Z'Z_(-1)) 

A_new = A;
A_new(1:r,:) = EZZ_FB(1:r,:)/(EZZ_BB);
Q_new = Q;
Q_new(1:r,1:r) = (EZZ(1:r,1:r) - A_new(1:r,:)*EZZ_FB(1:r,:)') / T;


%E(Y'Y) & E(Y'Z) 
nanY = isnan(y);
y(nanY) = 0;

denom = zeros(n*r,n*r);
nom = zeros(n,r);
for t=1:T
    nanYt = diag(~nanY(:,t));
    denom = denom + kron(Zsmooth(1:r,t+1)*Zsmooth(1:r,t+1)'+Vsmooth(1:r,1:r,t+1),nanYt);
    nom = nom + y(:,t)*Zsmooth(1:r,t+1)';
    
end

vec_C = (denom)\nom(:);
C_new = reshape(vec_C,n,r);


R_new = zeros(n,n);
for t=1:T
    nanYt = diag(~nanY(:,t));
    R_new = R_new + (y(:,t)-nanYt*C_new*Zsmooth(1:r,t+1))*(y(:,t)-nanYt*C_new*Zsmooth(1:r,t+1))'...
        +nanYt*C_new*Vsmooth(1:r,1:r,t+1)*C_new'*nanYt...
        +(eye(n)-nanYt)*R*(eye(n)-nanYt);
end

R_new = R_new/T;
RR = diag(R_new); %RR(RR<1e-2) = 1e-2; 
R_new = diag(RR);


% Initial conditions
Z_0 = Zsmooth(:,1); %zeros(size(Zsmooth,1),1); %
V_0 = Vsmooth(:,:,1);
%--------------------------------------------------------------------------

function [converged, decrease] = em_converged(loglik, previous_loglik, threshold, check_increased)
% EM_CONVERGED Has EM converged?
% [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
%
% We have converged if the slope of the log-likelihood function falls below 'threshold', 
% i.e., |f(t) - f(t-1)| / avg < threshold,
% where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
% 'threshold' defaults to 1e-4.
%
% This stopping criterion is from Numerical Recipes in C p423
%
% If we are doing MAP estimation (using priors), the likelihood can decrase,
% even though the mode of the posterior is increasing.

if nargin < 3, threshold = 1e-4; end
if nargin < 4, check_increased = 1; end

converged = 0;
decrease = 0;

if check_increased
    if loglik - previous_loglik < -1e-3 % allow for a little imprecision
        fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
        decrease = 1;
    end
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end



