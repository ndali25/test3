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

initialcond=zeros(T-p,n,draws);

for k=1:draws % draws
    for t=1:T-p % time
        for i=1:n % variables 
            
            initialcond(t,i,k)=Y(t,i)-sum(histdec(t,:,i,k),2);  
            
        end   
    end
end