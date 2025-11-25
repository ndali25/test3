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
var_pos = [1,1,1,1,1];       % position of corresponding variable to shock 
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

        %     Check magnitude restriction
        if magnitude==1 && jj==1
            if (chk(3,1)-chk(1,1)) > 0 
                jj = 0;
                break
            end
        end
        if magnitude==1 && jj==2
            if (chk(3,1)-chk(1,1)) < 0
                jj = 0;
                break
            end
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



%%

Bdraw_all_cFirst = Bdraw_all;

BigA=zeros(n*p,n*p,draws);
Sigma=zeros(n,n,draws);
errornorm=zeros(T-p,n,draws);
fittednorm=zeros(T-p,n,draws);
for i=1:draws 
    [Y,X,Y_initial]     = SUR(data,p,c);
    BigA(:,:,i)         = [Bdraw_all_cFirst(1+c:end,:,i)'; eye(n*p-n) zeros(n*p-n,n)]; % (n*p)x(n*p) matrix
    if abs(eig(BigA(:,:,i)))>=1
        disp('Unstable')
        Bdraw_all_cFirst(:,:,i) = [];
        BigA(:,:,i)= [];
    end  
    errornorm(:,:,i)    = Y-X*Bdraw_all_cFirst(:,:,i);
    fittednorm(:,:,i)   = X*Bdraw_all_cFirst(:,:,i);
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

candidateirf_cum = cumsum(candidateirf,3);

%% Reshape the matrices into a 3D object:

% For each draw, compute a matrix with the IRFs for each variable and each
% shock for the entire horizon considered, i.e. hor x n*n for the # of
% draws:

candidateirf_wold=zeros(hor,n*n,draws); 

for k=1:draws
    
candidateirf_wold(:,:,k)=(reshape(permute(candidateirf(:,:,:,k),[3 2 1]),hor,n*n,[]));
candidateirf_wold_cum(:,:,k)=(reshape(permute(candidateirf_cum(:,:,:,k),[3 2 1]),hor,n*n,[]));

end