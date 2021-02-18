function [Xs,slice_rets] = solve_second_twosided(T,slices,b,tolerance)
% Iteratively solves the linear system in step 3 of the algorithm

s = size(T);
n = s(1);
maxiter = 5000;
printiter = 50;
stop_tol = tolerance;


slice_list = 1:2:slices;


binds = ((1:slices)-1)*(2*b-1)+1;


bad_rows = cell(slices,1);

slice_nums = binds;
slice_rets = binds(slice_list);

for num = slice_list
    bad_rows{num} =  max( binds(num) - 2*(b-1) ,1):min(binds(num)+(2*b-2),n);

end

As = cell(slices,1);
for num = slice_list
    As{num} = squeeze(T(:,:,slice_nums(num)));

end


Xs = cell(slices,1);

for num = slice_list
    Xs{num} = zeros(n,n);
end

LHS = cell(slices,slices);

for num = slice_list
    for nnum = slice_list
        LHS{num,nnum} = (As{nnum}'*As{num} - As{num}'*As{nnum});

    end
end


for iter = 1:maxiter
    Xsprev = Xs;

    for num = slice_list
        Xs{num} = solve_sub_slice_second_twosided(As,LHS,Xs,num,b,binds,bad_rows,slice_list);
    end
    

        
    stop = calc_stop(Xs,Xsprev);
    
    if mod(iter,printiter) == 0
        fprintf('iter: %d, rel. change %e, tol. %e \n', iter, stop, stop_tol)
    end
    if stop < stop_tol
        break
    end
end
