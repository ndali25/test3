function [f,sr,var_pos] = var_restrictions(identification)
    if strcmp(identification,'sign')==1 % Sign restrictions on impact
        % No zero restrictions
                    %supply,    demand
                    % Short term zero restrictions
        f      = [  1,          1; %GDP
                    1,          1; %Prices
                    % Long term zero restrictions
                    1,          1; %GDP
                    1,          1]; %Prices
        
        % Sign restrictions
                    %supply,    demand
        sr      = [ +1,         +1;  %GDP
                    -1,         +1];  %Prices

        var_pos = [1,1];       % position of corresponding variable to shock
    elseif strcmp(identification,'cholesky')==1 % Cholesky
        % Zero restriction
                        %supply,    demand
                        % Short term zero restrictions
        f          = [  0,          1; %GDP
                        1,          1; %Prices
                        % Long term zero restrictions
                        1,          1; %GDP
                        1,          1]; %Prices
        
        % Sign restrictions (only for normalization)
                    %supply,    demand
        sr      = [ nan,        nan;  %GDP
                    nan,        nan];  %Prices

        var_pos = [1,1];       % position of corresponding variable to shock

    elseif strcmp(identification,'cholesky_sim')==1 % Cholesky
        % Zero restriction
                        %shock1,    shock2
                        % Short term zero restrictions
        f          = [  1,          0; %Y1
                        1,          1; %Y2
                        % Long term zero restrictions
                        1,          1; %Y1
                        1,          1]; %Y2
        
        % Sign restrictions 
                    %supply,    demand
        sr      = [ nan,        nan;  %Y1
                    nan,         nan];  %Y2

        var_pos = [1,2];       % position of corresponding variable to shock
    
    elseif strcmp(identification,'bq')==1 % Blanchard-Quah
        % Zero restriction
                        %supply,    demand
                        % Short term zero restrictions
        f          = [  1,          1; %GDP
                        1,          1; %Prices
                        % Long term zero restrictions
                        1,          0; %GDP
                        1,          1]; %Prices
        
        % Sign restrictions 
                    %supply,    demand
        sr      = [ nan,        nan;  %GDP
                    nan,        nan];  %Prices

        var_pos = [1,1];       % position of corresponding variable to shock
    else
        error('identification should be sign, cholesky or bq')    
    end 

end



