function xFS = runKF_DK(y, A, C, Q, R, x_0, Sig_0,c1,c2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs Kalman filter and smoother

% run the filter
S = SKF(y,C,R,A,Q, x_0, Sig_0,c1,c2);
% run the smoother
S = FIS(y,C,R,A,S); 

xFS = S.AmT;


%______________________________________________________________________
function S = SKF(Y,Z,R,T,Q,A_0,P_0,c1,c2)
%______________________________________________________________________
% Kalman filter for stationary systems with time-varying system matrices
% and missing data.
%
% The model is        y_t   = Z * a_t + eps_t       
%                     a_t+1 = T * a_t + u_t       
%
%______________________________________________________________________
% INPUT  
%        Y         Data                                 (nobs x n)  
% OUTPUT 
%        S.Am       Predicted state vector  A_t|t-1      (nobs x m)  
%        S.AmU      Filtered  state vector  A_t|t        (nobs+1 x m)  
%        S.Pm       Predicted covariance of A_t|t-1      (nobs x m x m)  
%        S.PmU      Filtered  covariance of A_t|t        (nobs+1 x m x m)  
%        S.loglik   Value of likelihood function
  
% Output structure & dimensions
  [n, m] = size(Z);
  nobs  = size(Y,2);
  
  S.Am  = nan(m,nobs);   S.Pm  = nan(m,m,nobs);
  S.AmU = nan(m,nobs);   S.PmU = nan(m,m,nobs);
  S.ZF = cell(nobs);
  S.V = cell(nobs);
  %______________________________________________________________________
  Au = A_0;  % A_0|0;
  Pu = P_0;  % P_0|0
  
  for t = 1:nobs
%       t
      % A = A_t|t-1   & P = P_t|t-1

      A   = T*Au+c2;
      P   = T*Pu*T' + Q;
      P   =  0.5 * (P+P');
      
      % handling the missing data
      [y_t,Z_t,R_t,c1_t] = MissData(Y(:,t),Z,R,c1);

      if isempty(y_t)
          Au = A;
          Pu = P;
          ZF = zeros(m,0);
          V = zeros(0,1);

      else
          PZ = P*Z_t';
          F  = (Z_t*PZ + R_t);
          ZF = Z_t'/F;
          PZF = P*ZF;

          V   = y_t - Z_t*A-c1_t;
          Au  = A  + PZF * V;
          Pu  = P  - PZF * PZ';
          Pu   =  0.5 * (Pu+Pu');
      end
      S.ZF{t} = ZF;
      S.Am(:,t)   = A;
      S.Pm(:,:,t) = P;
      S.V{t} = V;
      
      S.AmU(:,t)    = Au;
      S.PmU(:,:,t)  = Pu;
  end % t
 

%______________________________________________________________________
function S = FIS(Y,Z,R,T,S)
%______________________________________________________________________
% Fixed interval smoother (see Durbin and Koopman, 2001, p. 64-71)
% FIS returns the smoothed state vector AmT and its covar matrix PmT             
% Use this in conjnuction with function SKF
%______________________________________________________________________
% INPUT  
%        Y         Data                                 (nobs x n)  
%        S Estimates from Kalman filter SKF                                                          
%          S.Am   : Estimates     a_t|t-1                  (nobs x m) 
%          S.Pm   : P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
% OUTPUT 
%        S Smoothed estimates added to above
%          S.AmT  : Estimates     a_t|T                    (nobs x m) 
%        where m is the dim of state vector and t = 1 ...T is time

  [m, nobs]        = size(S.Am);
  S.AmT           = zeros(m,nobs);
  S.PmT           = zeros(m,m,nobs);
  S.AmT(:,nobs)   = squeeze(S.AmU(:,nobs))  ;

  r = zeros(m,1);

  for t = nobs:-1:1
	[y_t,Z_t] = MissData(Y(:,t),Z,R,zeros(length(Y(:,t)),1));

	r = S.ZF{t}*S.V{t}+(T*(eye(m)-squeeze(S.Pm(:,:,t))*S.ZF{t}*Z_t))'*r;
	S.AmT(:,t) = S.Am(:,t)+ squeeze(S.Pm(:,:,t))*r;
      
     
  end

%______________________________________________________________________
function [y,C,R,c1]  = MissData(y,C,R,c1)
%______________________________________________________________________
% PROC missdata                                                        
% PURPOSE: eliminates the rows in y & matrices C, R and vector c1 that correspond to     
%          missing data (NaN) in y                                                                                  
%______________________________________________________________________
  ix = ~isnan(y);

  y  =  y(ix);
  c1 =  c1(ix);
  C  =  C(ix,:);  
  R  =  R(ix,ix);
