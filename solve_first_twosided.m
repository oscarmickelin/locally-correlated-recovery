function [Xs,slice_nums,iter] = solve_first_twosided(T,slices,b,tolerance)
% Iteratively solves the linear system in step 2 of the algorithm
s = size(T);
n = s(1);
maxiter = 5000;
printiter = 50;
stop_tol = tolerance;


binds = ((1:slices)-1)*(2*b-1)+1;

bad_rows = cell(slices,1);

slice_nums = binds;

for num = 1:slices
    bad_rows{num} =  max(((num-1)*(2*b-1) + 1 - 2*(b-1) ),1):num*(2*b-1);
end

As = cell(slices,1);
for num = 1:slices
    As{num} = squeeze(T(:,:,slice_nums(num)));
end

Xs = cell(slices,1);

for num = 1:slices
    Xs{num} = zeros(n,n);
end

LHS = cell(slices,slices);

for num = 1:slices
    for nnum = 1:slices
        LHS{num,nnum} = As{nnum}*As{num}' - As{num}*As{nnum}';

    end
end

for iter = 1:maxiter

    Xsprev = Xs;

    for num = 1:slices
        Xs{num} = solve_sub_slice_twosided(As,LHS,Xs,num,b,binds);
    end
        
    stop = calc_stop(Xs,Xsprev,bad_rows);

    if mod(iter,printiter) == 0
        fprintf('iter: %d, rel. change %e, tol. %e \n', iter, stop, stop_tol)
    end
    if stop < stop_tol
        break
    end
end
