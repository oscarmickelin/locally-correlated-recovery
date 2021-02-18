function mat = solve_sub_slice_second_twosided(As,LHS,Xs,num,b,binds,bad_rows,slice_list)
%Updates one matrix X_num in step 3 of the algorithm

mat = Xs{num};

[n,~] = size(As{1});

slices = length(As);
As_each = cell(slices-1,1);
LHS_each = cell(slices-1,1);
rows_avoid_each = cell(slices-1,1);


for nnum = slice_list
    
    if nnum <= num-1

        ind_nnum = bad_rows{nnum};

        rows_avoid_each{nnum} = ind_nnum;
        As_each{nnum} = -As{nnum}';
        LHS_each{nnum} = LHS{num,nnum} - (As{num}'*(Xs{nnum}) - Xs{nnum}'*(As{num}));

    elseif nnum >= num+1
        ind_nnum = bad_rows{nnum};

        rows_avoid_each{nnum-2} = ind_nnum;
        As_each{nnum-2} = -As{nnum}';
        LHS_each{nnum-2} = LHS{num,nnum} - (As{num}'*(Xs{nnum}) - Xs{nnum}'*(As{num})); 

    end
end



for k = 1:n
    mat(:,k) = solve_each_slice_second_twosided(As_each,LHS_each,mat,k,b,rows_avoid_each,binds,num,slice_list);
end

end

