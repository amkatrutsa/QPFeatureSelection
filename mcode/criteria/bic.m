function [ crit ] = bic( X, y, w, par )
% Compute BIC criteria for model. It computes through 
% BIC = RSS + p * log(m), where RSS is the residual sum of squares, 
% p is the number of the selected predictors, m is the number of objects. 
%
% Input:
% X - [m, p] - design matrix with shrinkage number of predictors
% y - [m, 1] - target vector
% w - [p, 1] - vector of parameters, getting from algorithm
% par - structure with additional parameters, in this caseis empty
%
% Output:
% crit - [1,1] - value for BIC criteria
%
% Author: Alexandr Katrutsa, 2016
% E-mail: aleksandr.katrutsa@phystech.edu
  
if(isempty(X))
    crit = Inf;
    return
end
rss = sumsqr(y - X * w);
crit = rss + size(X, 2) * log(size(X, 1));

end

