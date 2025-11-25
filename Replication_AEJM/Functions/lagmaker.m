function lags=lagmaker(y,p)

% Function to make a matrix of p lags of y:
% Author: Nicolo' Maffei Faccioli

    X=lagmatrix(y,1:p); % matrix of p lagged observations of y 
    X(1:p,:)=[]; % exclude the first p observations
    lags=X;

end