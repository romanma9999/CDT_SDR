function [mean_val, std_val] = calc_SDRT(S,N,nbits, experiments)
 
SDRT_sparsity = zeros(1,experiments);
for j = 1:experiments
     SDRT  = zeros(1,N);
     for i = 1:S         
        SDR{i} = zeros(1,N);
        SDR{i}(randperm(N,nbits)) = 1;
        SDRT = SDRT | SDR{i};
     end
     SDRT_sparsity(j) = sum(SDRT)/N;
 end
  
mean_val = mean(SDRT_sparsity)*100;
std_val =  std(SDRT_sparsity)*100;
end

