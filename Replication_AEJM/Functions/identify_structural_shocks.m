function R = identify_structural_shocks(Bdraw_all,Sigmadraw_all,f,sr,n,p,draws,noConstant,var_pos)
    if nargin<9
        var_pos = [1,1];       % position of corresponding variable to shock 
    end
    % Identify shocks
    % function identify_shocks
    len_IRF = 500; % length of impulse response function.  Set to 500 to check to see
    % long-run restrictions are satisfied.
    
    shocks = eye(n);
    
    [Q,index,flag] = findQs(n,f);
    
    if flag == 1
        error('Rank condition not satisfied, the model is overidentified');
    end
    k = n;
    shock_pos = logical(shocks); % position of shock
%     var_pos = [1,1];       % position of corresponding variable to shock 
    % eg monetary policy shock should result in a positive increase in interest 
    % rates, aggregate demand shock should result in a positive increase in gdp
    % etc.
    
    R = zeros(k,len_IRF,draws,length(shocks)); % Contains impulse resonse functions
    % 1st dimension = variable
    % 2nd dimension = time
    % 3rd dimension = draw
    % 4th dimension = shock
    
    counter = 1;
    counter_show = 0;
    
    
    sizemcmc=size(Bdraw_all,3);
     
        
    while counter < draws+1
        record = randi(sizemcmc);
        B = Bdraw_all(:,:,record);
    
        if noConstant ==1
            Btilde = B';
        else
            Btilde = B(2:end,:)';
        end
        Sigma = Sigmadraw_all(:,:,record);
    
        alpha = [Btilde;eye(k*(p-1)),zeros(k*(p-1),k)]; % Build companion form matrix
        C1 = chol(Sigma,'lower');
    
    
        C = generateDraw(C1,k);
        
        if noConstant ==1
            P = findP_noConst(C,B,Q,p,k,index);
            c = 0;
        else
            P = findP(C,B,Q,p,k,index);
            c = 1;
        end
        
        W = C*P;
    
        for jj = 1:length(shocks)
            
            if W(var_pos(jj),jj) < 0
                shock = -shocks(:,jj);
            else
                shock = shocks(:,jj);
            end
            
            V = zeros(k*p,len_IRF);
            
            V(1:k,1) = W*shock;
            
            chk = W*shock;
            sr_index = ~isnan(sr(:,jj));
            tmp = sign(chk(sr_index)) - sr(sr_index,jj);
    
            if any(tmp~=0)
                jj = 0;
                break
            end
            
            for ii = 2:len_IRF
                V(:,ii) = alpha*V(:,ii-1);
            end
            
            R(:,:,counter,jj) = V(1:k,:);
            
        end
        
        if jj == length(shocks)
            counter = counter + 1;
            counter_show = counter_show +1;
            if counter_show==100
                disp(counter-1)
                counter_show=0;
            end
        end
        
    end
end