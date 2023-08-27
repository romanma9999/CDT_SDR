clc

clear all
close all

CDT_type = 1; % 0 is substructive, 1 is additive, 2 is combined
CDT_string{1} = "Substructive CDT: ";
CDT_string{2} = "Additive CDT: ";
CDT_string{3} = "Stable CDT: ";

experiments = 50;
N = 2000;
nbits = 40;

sparsity = floor(100*nbits/N)/100;

%SS = [2:16];
SS = [2:10];
NS = length(SS);
expected_residual_ratio = zeros(1,NS);
residual_sdr_mean_vec = zeros(1,NS);
residual_sdr_std_vec = zeros(1,NS);
final_sdr_mean_vec = zeros(1,NS);
final_sdr_std_vec = zeros(1,NS);
sdr_mean_vec = zeros(1,NS);
sdr_std_vec = zeros(1,NS);
k_mean_vec = zeros(1,NS);
k0_mean_vec = zeros(1,NS);
k1_mean_vec = zeros(1,NS);
k_std_vec = zeros(1,NS);
k_vec = zeros(1,NS);
for is = 1:length(SS)
    PKZT = zeros(1,N);
    S = SS(is);
    disp(['calc for S=' num2str(S)]);
    [sdr_mean_vec(is), sdr_std_vec(is)] = calc_SDRT(S,N,nbits,experiments);
    SDRT = zeros(1,N);
    while(1)
        SDR_tmp = zeros(1,N);
        for i = 1:S
            SDR{i} = zeros(1,N);
            SDR{i}(randperm(N,nbits)) = 1;
            SDR_tmp = SDR_tmp | SDR{i};
        end
        SDRT_sparsity = 100*sum(SDR_tmp)/N;
        if (abs(SDRT_sparsity - sdr_mean_vec(is)) < 0.2)
            SDRT = SDR_tmp;
            break;
        end
    end    

    [SDRF,NKF0,NKF1] = CDT(experiments,sparsity,N,SDRT,CDT_type);
     
    rs = zeros(1,experiments);
    fs = zeros(1,experiments);
    expected_residual_ratio(is) = nbits / S;
    for j = 1:experiments
      fs(j) = sum(SDRF{j})/N;
      for idx_sdr = 1:S
        rs(j) = rs(j) + sum(SDRF{j} & SDR{idx_sdr});
      end
      rs(j) = rs(j)/(S);
     end
     residual_sdr_mean_vec(is) = mean(rs);   
     residual_sdr_std_vec(is) = std(rs);   
     final_sdr_mean_vec(is) = mean(fs).*100;
     final_sdr_std_vec(is) = std(fs).*100;
    
    k_mean_vec(is) = mean(cell2mat(NKF0) + cell2mat(NKF1));
    k_std_vec(is) = std(cell2mat(NKF0) + cell2mat(NKF1));
    k0_mean_vec(is) = mean(cell2mat(NKF0));
    k1_mean_vec(is) = mean(cell2mat(NKF1));
     if CDT_type == 0
        k_vec(is) = log(S)/(sparsity*S);
    end
    if CDT_type == 1
        k_vec(is) = log(1.0 - 1.0/S)/log(1.0 - sparsity*S);
    end
end

figure;
plot_mean_var(sdr_mean_vec,sdr_std_vec);
hold on;
plot(SS, SS.*sparsity*100,'--r','LineWidth', 2);
title(['Combined SDR sparsity before CDT, for s = '  num2str(sparsity*100) '% and N = ' num2str(N)]);
xlabel('#Combined SDRs');
ylabel('Sparsity (%)');
legend("Combined SDR sparsity std","Combined SDR sparsity mean","sparsity * #SDRs");

figure;
plot_mean_var(k_mean_vec,k_std_vec);

if CDT_type == 0
    hold on;
    plot(SS, k_vec,'--r','LineWidth', 2);
    legend("K std","K mean","K expected");
    title(CDT_string{CDT_type+1} +'Number (K) of Permutations as function of combined #SDR to achieve s = ' +num2str(sparsity*100) +'%, ' +num2str(experiments) +' runs');
end
if CDT_type == 1
    hold on;
    plot(SS, k_vec,'--r','LineWidth', 2);
    legend("K std","K mean","K expected");
    title(CDT_string{CDT_type+1} + ' Number (K) of Permutations as function of combined #SDR to achieve s = ' + num2str(sparsity*100) +' %, '+ num2str(experiments) +' runs');
end
if CDT_type == 2
    hold on;
    plot(SS,k0_mean_vec,':','LineWidth', 2);
    hold on;
    plot(SS,k1_mean_vec,'-.','LineWidth', 2);
    legend("K total std","K total mean","K substructive","K additive");
    title(CDT_string{CDT_type+1} + ' Number (K) of Permutations as function of combined #SDR to achieve s = ' + num2str(sparsity*100) + '%, ' + num2str(experiments) + ' runs');
end

figure;
plot_mean_var(final_sdr_mean_vec,final_sdr_std_vec);
hold on;
title(CDT_string{CDT_type+1} + ' sparsity after CDT, s = ' + num2str(sparsity*100) + '% and N = ' + num2str(N));
xlabel('#Combined SDRs');
ylabel('Sparsity (%)');
legend("Final SDR sparsity std","Final SDR sparsity mean");

figure;
plot_mean_var(residual_sdr_mean_vec,residual_sdr_std_vec);
hold on;
plot(SS, expected_residual_ratio,'--r','LineWidth', 2);
title(CDT_string{CDT_type+1} + 'Average ON bits taken from each SDR when using CDT, s = '  + num2str(sparsity*100) + '% and N = ' + num2str(N));
xlabel('#Combined SDRs');
ylabel('#SDR ON bits');
legend("SDR ON bits std","SDR ON bits mean","expected ratio");


