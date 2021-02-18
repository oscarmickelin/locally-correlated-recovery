function [Asolve, Bsolve] = jennrich(S,r)
%%% Implements Jennrich's algorithm for orthogonally decomposable tensors
s = size(S);
n = s(3);

z1 = randn(n,1);
z1 = z1/norm(z1);

M1 = 0;

for i = 1:n
    M1 = M1 + S(:,:,i)*z1(i);
end

[Asolve, sings, Bsolve] = svd(M1, 'econ');


end
