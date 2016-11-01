function [ w ] = ElasticNet(X, y, par)
% Function wraps the elastic net regularization of the OLS
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

[W, inf] = lasso(X, y, 'Alpha', 0.5, 'CV', 5);
fprintf('Corresponding lambda = %d\n', inf.LambdaMinMSE);
w = W(:, inf.IndexMinMSE);
end

