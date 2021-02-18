%%% Example that bounds the possible values of b as function of the number
%%% of slices used
clc; clf; clear all; close all

n = 100;
bs = 1:8;
slices = 5:1:9;
tolerance = 1e-6;
errs = zeros(length(bs)-1,max(slices));
iters = zeros(length(bs)-1,max(slices));

repeats = 1;

for n = n
    r = n;
    for b = bs
        disp(b)
        by = 1;
        bx = b;
        for slice = slices
            if (2*b + (slice-2)*(2*b-1)) > n
                err = 1;
            else
                [err,iter] = solve_unknown_twosided(n,r,b,slice,tolerance);
            end
            errs(b,slice-1) = err;
        end
    end
end

errs




%%% plot
toPlt = errs(:, 4:8) < 0.03;
toPlt = [toPlt; zeros(1,5)];
imagesc(flipud(toPlt'))

ln = 2;

line(repmat([0.5:1:9.5],2,1),repmat([0;7], 1, 10), 'Color', [0 0 0], 'LineWidth', ln); 
line(repmat([0;10], 1, 7), repmat([0.5:1:6.5],2,1), 'Color', [0 0 0], 'LineWidth', ln); 
newmap = contrast(toPlt);
axis equal

col1 = [0.4660, 0.6740, 0.1880];
col2 = [0.6350, 0.0780, 0.1840];
map = [col2; col1];

colormap(newmap)
colormap(map)
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20; 
xlabel('$$b$$', 'interpreter', 'latex', 'FontSize', 32)
ylabel('Number of slices', 'interpreter', 'latex', 'FontSize', 32)
yticklabels({'9','8','7','6','5'})
yticks(1:6)
xticks(1:9)
xlim([0.5,9.5])
ylim([0.5,5.5])

xticklabels({'0','1','2','3','4','5','6','7','8'})


%hack to get legend entries
hold on;
h = zeros(2, 1);

h(1) = plot(NaN,NaN, 'sk', 'MarkerSize', 20, 'MarkerFaceColor', col1, 'LineWidth', ln);
h(2) = plot(NaN,NaN, 'sk', 'MarkerSize', 20, 'MarkerFaceColor', col2, 'LineWidth', ln);

legend(h, 'Recovered','Not recovered', 'interpreter', ...
    'latex', 'FontSize', 32, 'location', 'NorthOutside', 'numcolumns', 2);
legend boxoff


pbaspect([2,1,1])
print -dpdf plot_bound_bs.pdf