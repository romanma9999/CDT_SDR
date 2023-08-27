function [outputArg1,outputArg2,outputArg3] = CDT(experiments,target_sparsity,N,SDRT,CDT_type)
 tmp_CDT_type = CDT_type;
 Peye = eye(N);
 Pk = Peye(randperm(N),:);
 PKZ = SDRT;
 for i = 1:experiments
      if tmp_CDT_type == 1 || tmp_CDT_type == 2
           SDR_FINAL = zeros(1,N);
      else
          SDR_FINAL = SDRT;
      end
    NK0 = 0;
    NK1 = 0;
    current_CDT_type = tmp_CDT_type;
     while (1)         
          Pk = Peye(randperm(N),:);
          PKZ = PKZ*Pk;
           if current_CDT_type == 1 || current_CDT_type == 2
                SDR_FINAL = SDR_FINAL | SDRT&PKZ;
                NK1 = NK1 + 1;
                current_sparsity = sum(SDR_FINAL)/N;
                if current_CDT_type == 1  && current_sparsity > target_sparsity 
                     break;
                end
                if current_CDT_type == 2 && current_sparsity > target_sparsity 
                     current_CDT_type = 0;
                end
           else
                SDR_FINAL = SDR_FINAL&(~PKZ);
                NK0 = NK0 + 1;
                 if sum(SDR_FINAL)/N < target_sparsity
                        break
                end
           end        
     end
     SDRF{i} = SDR_FINAL;
     NKF0{i} = NK0;
     NKF1{i} = NK1;
 end
outputArg1 = SDRF;
outputArg2 = NKF0;
outputArg3 = NKF1;
end

