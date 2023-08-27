function plot_mean_var(mean_vec,stddec_vec)

y = mean_vec; % your mean vector;
x = 2:numel(y)+1;
curve1 = y + stddec_vec;
curve2 = y - stddec_vec;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

fill(x2, inBetween, [.7 .7 .7]);
hold on;
plot(x, y, 'b', 'LineWidth', 2);
grid on

end

