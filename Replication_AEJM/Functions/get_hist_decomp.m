function [histdec, det_comp] = get_hist_decomp(data,p,c,Bdraw_all,R)
    %%
    [T,n] = size(data);
    draws = size(Bdraw_all,3);
    hor = size(R,2);

    BigA=zeros(n*p,n*p,draws);
    
    errornorm=zeros(T-p,n,draws);
    fittednorm=zeros(T-p,n,draws);
    for i=1:draws 
        [Y,X,~]     = SUR(data,p,c);
        BigA(:,:,i)         = [Bdraw_all(1+c:end,:,i)'; eye(n*p-n) zeros(n*p-n,n)]; % (n*p)x(n*p) matrix
        if abs(eig(BigA(:,:,i)))>=1
            disp('Unstable')
            Bdraw_all(:,:,i) = [];
            BigA(:,:,i)= [];
        end  
        errornorm(:,:,i)    = Y-X*Bdraw_all(:,:,i);
        fittednorm(:,:,i)   = X*Bdraw_all(:,:,i);
    end
    
    %% Restructure IRF object
    candidateirf = nan(n,n,hor,draws);
    
    for h=1:hor
        
       for k=1:draws
           
           for variable=1:n
               
               for shock=1:n
                    candidateirf(variable,shock,h,k) = R(variable,h,k,shock);
               end
           end
       end
    end
    
    % candidateirf_cum = cumsum(candidateirf,3);
    
    
    
    %% Historical decomposition:
    
    Y=data(p+1:end,:); 
    err=errornorm; % reduced form residuals
    xi=zeros(size(errornorm,1),size(errornorm,2),size(errornorm,3));
    HH=zeros(T-p,n,n,draws); % store the impact IRFs ^t 
    histdec=zeros(T-p,n,n,draws);
    
    for k=1:draws % draw k
        
        xi(:,:,k)=(candidateirf(:,:,1,k)\err(:,:,k)')'; 
        
        for t=1:T-p
                
                BigC=(BigA(:,:,k))^(t-1);
                BigH=BigC(1:n,1:n)*candidateirf(:,:,1,k); % C * S * H
                
                for j=1:n % variable j
                    
                HH(t,:,j,k)=BigH(j,1:n);
                
                end
                
        end
        
        for t=1:T-p
           
            for i=1:n % shock i 
                
                for j=1:n % variable j
                
                  histdec(t,i,j,k)=HH(1:t,i,j,k)'*flipud(xi(1:t,i,k));
                  
                end
                
            end
            
        end
    
    end
    
    % Deterministic component:
    
    det_comp=zeros(T-p,n,draws);
    
    for k=1:draws % draws
        for t=1:T-p % time
            for i=1:n % variables 
                
                det_comp(t,i,k)=Y(t,i)-sum(histdec(t,:,i,k),2);  
                
            end   
        end
    end

end