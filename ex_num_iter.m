%%% Example that bounds the number of iterations required for convergence
%%% as function of b and r

clc; clf; clear all; close all

ns = 100;
tolerance = 1e-7;
bs = 1:8;

rlist = 70:10:100;
rl = length(rlist);

errs = zeros(length(bs),rl);
iters = zeros(length(bs),rl);
slices = 7;

repeats = 10;

rind = 1;

for r = rlist
    for n = ns
        for b = bs
            disp([r,b])
            by = 1;
            bx = b;
            for i = 1:repeats
                [err,iter] = solve_unknown_twosided(n,r,b,slices,tolerance);
                errs(b,rind) = errs(b,rind) + err/repeats;
                iters(b,rind) = iters(b,rind) + iter/repeats;
            end

        end
    end
    rind = rind + 1;
end


%%% Plot

markers = {'o-', 's--', 'x:', 'd-.'};
lines = {'-', '--', ':', '-.'};

for r = 1:length(rlist)
    semilogy(bs-1, iters(:,r), markers{r}, 'LineWidth', 2)
    hold on
end

        
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
xlabel('$$b$$', 'interpreter', 'latex', 'FontSize', 32)
ylabel('Iterations', 'interpreter', 'latex', 'FontSize', 32)

legend({'$$r = 70$$','$$r = 80$$','$$r = 90$$','$$r = 100$$'}, ...
    'interpreter', 'latex', 'FontSize', 20, 'location', 'northwest')
grid on
ylim([6,2000])
yticks([1,100,1000])
ax.YMinorGrid = 'on';
ax.YAxis.MinorTickValues = [1:10:100, 100:100:1000];

pbaspect([3,1,1])

legend boxoff 
print -dpdf plot_num_iter.pdf