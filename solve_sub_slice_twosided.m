function mat = solve_sub_slice_twosided(As,LHS,Xs,num,b,binds)
%Updates one matrix X_num in step 2 of the algorithm

mat = Xs{num};

[n,~] = size(As{1});

slices = length(As);
As_each = cell(slices-1,1);
LHS_each = cell(slices-1,1);
inds_avoid_each = cell(slices-1,1);

ind_num = max(((num-1)*(2*b-1) + 1 - 2*(b-1) ),1):min(num*(2*b-1),n);
for nnum = 1:num-1

    ind_nnum = max(((nnum-1)*(2*b-1) + 1 - 2*(b-1) ),1):min(nnum*(2*b-1),n);

    inds_avoid_each{nnum} = union(ind_num, ind_nnum);
    As_each{nnum} = -As{nnum};
    LHS_each{nnum} = LHS{num,nnum} - (As{num}*(Xs{nnum}') - Xs{nnum}*(As{num}')); 
end

for nnum = num+1:slices

    ind_nnum = max(((nnum-1)*(2*b-1) + 1 - 2*(b-1) ),1):min(nnum*(2*b-1),n);

    inds_avoid_each{nnum-1} = union(ind_num, ind_nnum);
    As_each{nnum-1} = -As{nnum};
    LHS_each{nnum-1} = LHS{num,nnum} - (As{num}*(Xs{nnum}') - Xs{nnum}*(As{num}')); 
end


for k = 1:n
    mat(k,:) = solve_each_slice_twosided(As_each,LHS_each,mat,k,b,inds_avoid_each,binds,num);
end

end
