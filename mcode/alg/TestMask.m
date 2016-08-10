function [error] = TestMask(mask, features, target, param)
% Function TestMask estimates parameter with OLS and computes computes 
% the error function used features from mask vector.
%
% Input:
% mask - features mask (column) or masks (matrix)
% features - [n, m] - design matrix, where rows are objects, and colums are features
% target - [m, 1] - target vector
% 
% Output:
% error - scalar or [1, size(mask,1)] - scalar (for row mask) and vector (for matrix mask)
% with errors in parameter vector estimation

error = size(mask, 1);
i = 1;
for maskCol = mask
    % error(i) = ResidualNorm(features(:, find(maskCol)), target, param);
    error(i) = SumCriteria(features(:, find(maskCol)), target, param);
    i = i + 1;
end
end

function [error] = ResidualNorm(X, y, param)
w = lscov(X, y);    
error = norm(y - X * w);
end

function [error] = SumCriteria(X, y, param)
error = 0;
w = lscov(X, y);
for l = 1:length(param.crit)
    error = error + feval(param.crit{l}, X, y, w, param);
end
end
