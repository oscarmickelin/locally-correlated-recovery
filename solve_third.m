function Csolve = solve_third(T,Asolve,Bsolve,bx,by,b)
%Solves step 5 in the algorithm

[n,r] = size(Asolve);
Csolve = zeros(n,r);

for num = 1:n
    
    minind = max(1,num-bx+1);
    maxind = min(n,num+by-1);

    
    bvec = vec(double(T(:,:,num)));
    mat = ones(n,n);
    
    mat(minind:maxind,:) = zeros(maxind-minind+1,n);
    mat(:,minind:maxind) = zeros(n,maxind-minind+1);
    
    for i = 0:(by-1)
        mat(i+1:n+1:1+n*(n-i)) = zeros(n-i,1);
    end
    
    for i = 0:(bx-1)
        mat(i*n+1:n+1:end) = zeros(n-i,1);
    end

    ve = vec(mat);
    inds = (ve ~= 0);

    dg_ind = vec(eye(r,r)) > 0;

    big_mat = kron(Bsolve, Asolve);
    
    
    
    big_mat = big_mat(inds, dg_ind);
    
    bvec = bvec(inds);
    
    Csolve(num,:) = value(big_mat\bvec);
    
end
    
    
    
    
