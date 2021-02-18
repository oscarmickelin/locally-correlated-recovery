function [err,iter] = solve_unknown_noise_twosided(n,r,b,slices,tolerance,noise)
%%% Sets up a noisy version of a randomly generated recovery instance and
%%% solves it using the algorithm of the paper

A = orth(randn(n,r));
B = orth(randn(n,r));
C = randn(n,r);

T = double(ktensor({A,B,C}));

E = randn(n,n,n);
E = E/norm(tensor(E));
E = noise*E*norm(tensor(T));


T = T+E;

if (2*b + (slices-2)*(2*b-1)) > n
    disp('b too large...')
    disp('quitting')
    return
end
        

disp('......done')


for k = 1:n
    tmp = zeros(n,n);
    minind = max(1,k-b+1);
    maxind = min(n,k+b-1);
    
    tmp(minind:maxind,:) = randn(maxind-minind+1,n);
    tmp(:,minind:maxind) = randn(n,maxind-minind+1);
    
    for i = 0:(b-1)
        tmp(i+1:n+1:1+n*(n-i)) = randn(n-i,1);
    end
    
    for i = 0:(b-1)
        tmp(i*n+1:n+1:end) = randn(n-i,1);
    end
    
    T(:,:,k) = T(:,:,k) - tmp;

end

disp('...solving first')

[Xs,slice_nums,iter] = solve_first_twosided(T,slices,b,tolerance);


%%%

count = 1;
for num = slice_nums
    T(:,:,num) = T(:,:,num) + Xs{count};
    count = count+1;
end


disp('......done')

%%%
disp('...solving second')


[Xs,slice_rets] = solve_second_twosided(T,slices,b,tolerance);


count = 1;
for num = slice_rets
    T(:,:,num) = T(:,:,num) + Xs{count};
    count = count+2;
end

disp('......done')


%%%
Tslice = zeros(n,n,length(slice_rets));
count = 1;
for num = slice_rets
    Tslice(:,:,count) = T(:,:,num);
    count = count + 1;
end

disp('...finding ak, bk')
[Asolve,Bsolve] = jennrich(Tslice);
disp('......done')

%%%

disp('...solving third')

Csolve = solve_third(T,Asolve,Bsolve,b,b,b);

disp('......done')

%%%

err = norm(ktensor({Asolve, Bsolve, Csolve}) - ktensor({A,B,C}))/norm(ktensor({A,B,C}));