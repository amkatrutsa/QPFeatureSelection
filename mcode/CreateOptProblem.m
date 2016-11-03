function [Q, b] = CreateOptProblem(X, y, sim, rel)
% Function generates matrix Q and vector b 
% which represent feature similarities and feature relevances
%
% Input:
% X - [m, n] - design matrix
% y - [m, 1] - target vector
% sim - string - indicator of the way to compute feature similarities,
% support values are 'correl' and 'mi'
% rel - string - indicator of the way to compute feature significance,
% support values are 'correl', 'mi' and 'signif'
%
% Output:
% Q - [n ,n ] - matrix of features similarities
% b - [n, 1] - vector of feature relevances
%
% Author: Alexandr Katrutsa, 2016 
% E-mail: aleksandr.katrutsa@phystech.edu

if strcmp(sim, 'correl')
    Q = corrcoef(X);
else
    if strcmp(sim, 'mi')
        Q = zeros(size(X, 2));
        for i = 1:size(Q, 2)
            for j = i:size(Q, 2)
                Q(i, j) = information(X(:, i)', X(:, j)');
%                 Q(i, j) = mutualinfo(X(:, i), X(:, j));
            end
        end
        Q = Q + Q' - diag(diag(Q));
    end
    lambdas = eig(Q);
    min_lambda = min(lambdas);
    if min_lambda < 0
        Q = Q - min_lambda * eye(size(Q, 1));
    end
end
if strcmp(rel, 'correl')
    b = abs(corr(X, y));
    return
end
if strcmp(rel, 'mi')
    b = zeros(size(X, 2), 1);
    for i = 1:size(X, 2)
        b(i) = information(y', X(:, i)');
%         b(i) = mutualinfo(y, X(:, i));
    end
    return
end
if strcmp(rel, 'signif')
   lm = fitlm(X, y);
   p_val = lm.Coefficients.pValue(2:end);
   idx_zero_coeff = find(abs(lm.Coefficients.Estimate(2:end)) < 1e-7);
   nan_idx = isnan(p_val);
   p_val(nan_idx) = ones(sum(nan_idx), 1);
   b = 1 - p_val ./ sum(p_val); 
   b(idx_zero_coeff) = zeros(length(idx_zero_coeff), 1);
   return 
end
    
end