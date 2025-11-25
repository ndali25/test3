function lagged = lagmatrix(data, lags)
    % Custom implementation of lagmatrix
    % data: (T x n) matrix, where T is the number of time points and n is the number of variables.
    % lags: a vector of positive integers indicating the lag orders to include.
    
    [T, n] = size(data);
    maxLag = max(lags);
    lagged = NaN(T, n * length(lags)); % Initialize with NaNs for alignment.
    
    % Create lagged columns
    for i = 1:length(lags)
        lag = lags(i);
        lagged((lag + 1):end, (i - 1) * n + 1:i * n) = data(1:(end - lag), :);
    end
end
