function [ w ] = Ridge( X, y, par )
% Function wraps ridge regression
%
% Input:
% X - [m, n] - design matrix
% y - [m, 1] - target vector
%
% Output:
% w - [n, 1] - optimal parameter vector
%
% Author: Alexandr Katrutsa, 2016
% E-mail: aleksandr.katrutsa@phystech.edu

k = 0:0.0001:0.05; % See for limit value of k
W = ridge(y, X, k);
rss = zeros(1, size(W, 2));
for i = 1:size(W, 2)
    rss(i) = sumsqr(y - X * W(:, i));
end

[~, idx_min] = min(rss);
w = W(:, idx_min);
end