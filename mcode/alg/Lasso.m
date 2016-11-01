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

[W, fitinfo] = lasso(X, y, 'CV', 5);
fprintf('Corresponding lambda = %d\n', fitinfo.LambdaMinMSE);
w = W(:, fitinfo.IndexMinMSE);
end

