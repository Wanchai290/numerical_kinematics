function inv = clc_pinv(X)
% not very optimized
% in the sense we could've changed the computations
% in the main.m file to use SVD
% but there was no time left
% based on this result :
%   https://www.johndcook.com/blog/2018/05/05/svd/

[U, S, V] = svd(X);
S_size = size(S);
SS = zeros(S_size(2), S_size(1));
for i = 1:S_size(1)
    for j = 1:S_size(2)
        x = S(i, j);
        if x ~= 0
            SS(j, i) = 1/x;
        end
    end
end
inv = V * SS * transpose(U);
end