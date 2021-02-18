clc; clf; clear all; close all

n = 65;
tolerance = 1e-12;
b = 5;
rlist = 45:10:65;
rl = length(rlist);
noises = 10.^(-5:-1);
errs = zeros(length(noises),rl);
iters = zeros(length(noises),rl);
slices = 7;%6;

repeats = 100;

rind = 1;
for r = rlist
    currind = 1;
    for noise = noises
        for i = 1:repeats
            disp([r, noise, i])
            [err,iter] = solve_unknown_noise_twosided(n,r,b,slices,tolerance,noise);
            errs(currind,rind) = errs(currind,rind) + err/repeats;
            iters(currind,rind) = iters(currind,rind) + iter/repeats;
        end
        loglog(noises, errs, 'bs-')
        pause(0.1)
        currind = currind + 1;
    end
    rind = rind + 1;
end




%%% Plot
clf
figure(1)
markers = {'o-', 's--', 'x:', 'd-.'};

markers_jenn = {'^-', 'v--', '>:', '<-.'};



for r = 1:length(rlist)
    loglog(noises, errs(:,r), markers{r}, 'LineWidth', 2)
    hold on
end

for r = 1:length(rlist)
    loglog(noises, errs_jenn(:,r), markers_jenn{r}, 'LineWidth', 2)
    hold on
end

ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
xlabel('Relative perturbation', 'interpreter', 'latex', 'FontSize', 32)
ylabel('Relative error', 'interpreter', 'latex', 'FontSize', 32)

legend({"$$r = 35$$ ","$$r = 45$$ ","$$r = 55$$ ","$$r = 65$$ "}, ...
    'interpreter', 'latex', ...
   'location', 'SouthEast', 'FontSize', 20,'NumColumns',2)


pbaspect([3,1,1])
grid on

legend boxoff
print -dpdf plot_error_tol.pdf