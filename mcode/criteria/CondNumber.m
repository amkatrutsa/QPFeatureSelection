function [ kappa ] = CondNumber( X, y, w, par )
% Compute conditional number for shrinking design matrix X. 
%
% Input:
% X - [m, p] - design matrix with shrinkage number of predictors
% y - [m, 1] - target vector
% w - [p, 1] - vector of parameters, getting from algorithm
% par - structure with additional parameters, in this case is empty
%
% Output:
% kappa - [1,1] - conditional number for given design matrix
%
% Author: Alexandr Katrutsa, 2016 
% E-mail: aleksandr.katrutsa@phystech.edu

if(isempty(X))
    kappa = Inf;
    return
end
A = X' * X;
d = eig(A);
kappa = log(abs(max(d)) / abs(min(d)));
end