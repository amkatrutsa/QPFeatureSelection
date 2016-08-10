function [w] = Lasso(X, y, par)
% Function wraps the Lasso regularization of the OLS
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

[W, fitinfo] = lasso(X, y);
[~, idx_minMSE] = min(fitinfo.MSE);
w = W(:, idx_minMSE);
end

