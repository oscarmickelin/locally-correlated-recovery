function res = calc_stop(Xs,Xsprev,bad_rows)
%%% Calculates the relative improvement between current iterate Xs and the
%%% previous iterate Xsprev, for elements contained in bad_rows

slices = length(Xs);
[n,~] = size(Xs{1});


stop = 0;
stopn = 0;


for num = 1:slices
    if nargin == 3
        inds = setdiff(1:n, bad_rows{num});
    else
        inds = 1:n;
    end
    if numel(Xs{num}) ~= 0
        stop = stop + norm(Xs{num}(inds,:) - Xsprev{num}(inds,:), 'fro')^2;
        stopn = stopn + norm(Xsprev{num}(inds,:), 'fro')^2;
    end
end

res = sqrt(stop/stopn);

end