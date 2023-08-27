clc
clear all
close all

theta = 0.5;
ww = [16 32 64];
ww = 1;
nn = 1400:200:2048;
res = zeros(length(ww),length(nn));
for j = 1:length(ww);
    for i = 1:length(nn)
        n = nn(i)
        %w = ww(j)
        w = n*0.02;
        tmp = 0;
        for b = theta*w:w
            tmp = tmp + nchoosek( w , b ) * nchoosek(n - w , w - b);
        end
        res(j,i) = tmp/nchoosek(n,w);
    end
end

figure;
for i = 1:length(ww)
plot(nn,res(i,:));
hold on;
end
xlabel('SDR length');
ylabel('P false positive');
grid on;