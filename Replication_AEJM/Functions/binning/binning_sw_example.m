%% Copyright Andrew Binning 2013
% Please feel free to use and modify this code as you see if fit. If you
% use this code in any academic work, please cite 
% Andrew Binning, 2013.
% "Underidentified SVAR models: A framework for combining short and long-run restrictions with sign-restrictions,"
% Working Paper 2013/14, Norges Bank.

%% Load data and estimate reduced form VAR model
% data order:
% 'robs','dc','dinve','dy','labobs','pinfobs','dw'
% interest rates, consumption, investment, gdp, hours, inflation, wages

clear all

[num,txt,raw] = xlsread('Smets_Wouters_data.xlsx');

data = num(:,2:end);

data = data(74:end,[1,4:end]);

p = 2; %set lags

Y = data(p+1:end,:);
X = xlags(data,p);

[n,k] = size(Y);

X = [ones(n,1),X];

B = (X'*X)\X'*Y;

u = Y - X*B;

Sigma = u'*u/(n-k);

%% Short run, long run and sign restrictions

% f = 2.k x k matrix contain ones and zeros. The zeros occur where the zero
%   restrictions are on the short and long-run impact matrices.  The top k x
%   k block are the short run restrictions, while the bottom k x k
%   restrictions are the long run restrictions.
% sr = k x k matrix of sign restrictions to impose on impact.
% mp = monetary policy shock
% ad = aggregate demand shock
% as = aggregate supply shock
% ls = labour markup shock

tic;

    %mp,ad,as,ls
f = [1,1,1,1,1; %interest rates
     1,1,1,1,1; %gdp
     1,1,1,1,1; %hours
     1,1,1,1,1; %inflation
     1,1,1,1,1; %wages
     %=========
     1,1,1,1,1; %interest rates
     0,0,1,0,0; %gdp
     1,1,1,1,1; %hours
     1,1,1,1,1; %inflation
     1,1,1,1,1];%wages


      %mp %ad %as %ls
sr = [ +1, +1, -1, -1,nan;  %interest rates
       -1, +1, +1, +1,nan;  %gdp
      nan,nan,nan, +1,nan;  %hours
       -1, +1, -1, -1,nan;  %inflation
      nan,nan,nan, -1,nan]; %wages

C1 = chol(Sigma,'lower');

draws = 1000; % Number of draws

T = 500; % length of impulse response function.  Set to 500 to check to see
% long-run restrictions are satisfied.

shocks = eye(k);

[Q,index,flag] = findQs(k,f);

if flag == 1
    error('Rank condition not satisfied, the model is overidentified');
end

shock_pos = logical(shocks); % position of shock
var_pos = [1,2,2,3,3];       % position of corresponding variable to shock 
% eg monetary policy shock should result in a positive increase in interest 
% rates, aggregate demand shock should result in a positive increase in gdp
% etc.

R = zeros(k,T,draws,length(shocks)); % Contains impulse resonse functions
% 1st dimension = variable
% 2nd dimension = time
% 3rd dimension = draw
% 4th dimension = shock

counter = 1;

Btilde = B(2:end,:)';
alpha = [Btilde;eye(k*(p-1)),zeros(k*(p-1),k)]; % Build companion form matrix

while counter < draws+1
    
    C = generateDraw(C1,k);
    
    P = findP(C,B,Q,p,k,index);
    
    W = C*P;

    for jj = 1:length(shocks)
        
        if W(var_pos(jj),jj) < 0
            shock = -shocks(:,jj);
        else
            shock = shocks(:,jj);
        end
        
        V = zeros(k*p,T);
        
        V(1:k,1) = W*shock;
        
        chk = W*shock;
        sr_index = ~isnan(sr(:,jj));
        tmp = sign(chk(sr_index)) - sr(sr_index,jj);

        if any(tmp~=0)
            jj = 0;
            break
        end
        
        for ii = 2:T
            V(:,ii) = alpha*V(:,ii-1);
        end
        
        R(:,:,counter,jj) = V(1:k,:);
        
    end
    
    if jj == length(shocks)
        counter = counter + 1
    end
    
end

toc;

plot(squeeze(R(1,1:20,:,1))); % Plot first 20 periods of interest rate shock after a monetary policy shock
