function [Y,X,Y_initial]=SUR_sur(data,p,c,y0_bar,delta)
% Function that computes Y and X of the SUR representation 
% Author: Nicolo' Maffei Faccioli
%
% Y = X*PI + e
%
% INPUTS: 
% Data: T x n dataset 
% p: number of lags 
% 
% OUTPUTS:
% Y: T x n matrix of the SUR representation
% X: T x (n*p) matrix of the SUR representation

Y=(data(p+1:end,:));

if c==1
    X =[ones(length(lagmatrix(data,1:p)),1) lagmatrix(data,1:p)]; 
else 
    X=lagmatrix(data,1:p);
end 
    
X(1:p,:)=[];   

% Discarded observations:

Y_initial=data(1:p,:);

Y_d = y0_bar./delta;
X_d = [1/delta, Y_d];


Y = [Y;Y_d];
X = [X;X_d];

end
