function v = solve_each_slice_second_twosided(As,LHS,X,k,b,rows_avoid_each,binds,curr_num,slice_list)
%%% solves one iteration of A*X - X'*A', where X is a matrix with band structure


[n,~] = size(As{1});
samp = length(As);


minind_dbl = max(1,binds(curr_num)-2*(b-1));
maxind_dbl = min(n,binds(curr_num)+2*(b-1));


inds_col = minind_dbl:maxind_dbl;

if mod(samp,2) == 0
    coeff = zeros((samp)/2*(n-1),length(inds_col));
    bvec = zeros((samp)/2*(n-1),1);
    maxind = samp;
else
    coeff = zeros((samp+1)/2*(n-1),length(inds_col));
    bvec = zeros((samp+1)/2*(n-1),1);
    maxind = samp-1;
end

offset = 1;

tmp = X;
tmp(:,k) = 0;


for nnum = 1:2:maxind

        inds = [];
    
        inds = setdiff(1:n,union(inds,k));
        
        tmp_inds = (offset):(offset+length(inds)-1);

        coeff(tmp_inds,:) = As{nnum}(inds,inds_col);  
        pert_LHS = LHS{nnum}(:,k) + tmp'*(As{nnum}(k,:)');

        bvec(tmp_inds) = pert_LHS(inds);
        offset = offset + length(inds);

end


if offset == 1
    v = X(:,k);
else
    coeff = coeff(1:offset-1,:);
    bvec = bvec(1:offset-1);
    vres = coeff\bvec;
    v = zeros(n,1);
    v(inds_col) = vres;


end


end
