function [W, H] = NMF(X, K, MAXITER)
%Euclidean distance
[m,n]= size(X);
% rand('seed',0)
W = rand(m, K);
H = rand(K, n);
for i=1:MAXITER
    H = H .* (W'*X)./(W'*W*H+eps) ;
    W = W .* (X*H')./(W*H*H'+eps);
end
end