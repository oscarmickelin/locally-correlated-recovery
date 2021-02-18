function v = solve_each_slice_twosided(As,LHS,X,k,b,inds_avoid_each,binds,curr_num)
%%% solves one iteration of A*X' - X*A', where X is a matrix with structure
%%% band-diagonal + cross

[n,~] = size(As{1});
samp = length(As);

bandinds = max(1,(binds(curr_num)-(b-1))):min(binds(curr_num)+b-1,n);
bcurr = length(bandinds);

[n,~] = size(As{1});

minind = max(1,k-(b-1));
maxind = min(n,k+(b-1));

minind_dbl = max(1,k-2*(b-1));
maxind_dbl = min(n,k+2*(b-1));

coeff = zeros(samp*(n-1),bcurr+maxind-minind+1);

bvec = zeros(samp*(n-1),1);
offset = 1;

tmp = X;
tmp(k,:) = 0;

for num = 1:samp

    inds = inds_avoid_each{num};
    if ismember(k, inds)

    else
        inds = union(inds, minind_dbl:maxind_dbl);
        inds = setdiff(1:n, union(k,inds));
        
        tmp_inds = (offset):(offset+length(inds)-1);
        coeff(tmp_inds,:) = As{num}(inds,[bandinds, minind:maxind]);  
                
        pert_LHS = LHS{num}(:,k) + tmp*(As{num}(k,:)');

        bvec(tmp_inds) = pert_LHS(inds);
        offset = offset + length(inds);
    end
end



if offset == 1
    v = X(k,:);
else
    coeff = coeff(1:offset-1,:);
    bvec = bvec(1:offset-1);

    vres = coeff\bvec;
    v = zeros(1,n);
    v(bandinds) = vres(1:bcurr);

    v(minind:maxind) = vres(bcurr+1:end);

end

end