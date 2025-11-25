%% Fry and Pagan:

% Fry and Pagan IRFs:

resp_median=prctile(candidateirf,50,4);

distance_med=zeros(n,n,hor,draws);

for h=1:hor
    
   for k=1:draws
       
       for variable=1:n
           
           for shock=1:n

                   distance_med(variable,shock,h,k)=((candidateirf(variable,shock,h,k)-resp_median(variable,shock,h))/nanstd(squeeze(candidateirf(variable,shock,h,:)))).^2;
                   
           end
           
       end
   end
   
end

distance=squeeze(nansum(sum(nansum(distance_med(1:n,1:n,:,:),1),2),3));

for k=1:draws
    
if distance(k)==0
    distance(k)=NaN;
end

end

n_selected_draws = 10;

[distance_min,solutions]=mink(distance,n_selected_draws);


solution=solutions(1);

candidateirf_wold_fp=reshape(permute(candidateirf(:,:,:,solution),[3 2 1]),hor,n*n,[]); % Median-target IRFs

solution=solutions(2);

candidateirf_wold_fp_2=reshape(permute(candidateirf(:,:,:,solution),[3 2 1]),hor,n*n,[]); % Median-target IRFs

solution=solutions(3);

candidateirf_wold_fp_3=reshape(permute(candidateirf(:,:,:,solution),[3 2 1]),hor,n*n,[]); % Median-target IRFs

solution=solutions(4);

candidateirf_wold_fp_4=reshape(permute(candidateirf(:,:,:,solution),[3 2 1]),hor,n*n,[]); % Median-target IRFs

solution=solutions(5);

candidateirf_wold_fp_5=reshape(permute(candidateirf(:,:,:,solution),[3 2 1]),hor,n*n,[]); % Median-target IRFs

solution=solutions(6);

candidateirf_wold_fp_6=reshape(permute(candidateirf(:,:,:,solution),[3 2 1]),hor,n*n,[]); % Median-target IRFs

solution=solutions(7);

candidateirf_wold_fp_7=reshape(permute(candidateirf(:,:,:,solution),[3 2 1]),hor,n*n,[]); % Median-target IRFs

solution=solutions(8);

candidateirf_wold_fp_8=reshape(permute(candidateirf(:,:,:,solution),[3 2 1]),hor,n*n,[]); % Median-target IRFs

solution=solutions(9);

candidateirf_wold_fp_9=reshape(permute(candidateirf(:,:,:,solution),[3 2 1]),hor,n*n,[]); % Median-target IRFs

solution=solutions(10);

candidateirf_wold_fp_10=reshape(permute(candidateirf(:,:,:,solution),[3 2 1]),hor,n*n,[]); % Median-target IRFs