function [ pvalue ] = ftest( X, y, w, par)
% Function tests the null hypothesis that 
% there isn't appropriate predictors in the matrix X 
% Input:
% X - [m, p] - matrix with data set
% y - [m, 1] - target vector
% w - [p, 1] - vector of parameters given by algorithm
% par - structure with additional parameters, in this case is empty
%
% Output:
% pvalue - [1, 1] - p-value given by F-test
%
% Author: Alexandr Katrutsa, 2016
% E-mail: aleksandr.katrutsa@phystech.edu

if(isempty(X))
    pvalue = Inf;
    return
end
rss = sumsqr(y - X * w);
p = sum(w ~= 0);
tss = sumsqr(mean(y) - y);
if (rss > tss)
    pvalue = inf;
    return
end    
n = size(X, 1);
fstat = (abs(tss - rss) / p) / (rss / (n - p - 1));
pvalue = 2 * min(fcdf(fstat, p, n - p - 1), fcdf(1 / fstat, n - p - 1, p));
end

