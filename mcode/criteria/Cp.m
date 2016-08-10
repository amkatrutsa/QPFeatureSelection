function [ res ] = Cp( X, y, w, par)
% Compute Mallow's Cp criteria for model. 
% It computes through Cp = RSS_p / RSS_all + 2*p - m, 
% where RSS_p is the residual sum of squares for p selected predictors, 
% RSS_all is the residual sum of squares for all predictors, 
% p is the number of the selected predictors, m is the number of objects. 
% 
% Input:
% X - [m, p] - design matrix with shrinkage number of predictors
% y - [m, 1] - target vector
% w - [p, 1] - vector of parameters, getting from algorithm
% par - structure with additional parameters:
%           par.rss[1,1] --- RSS_all
% Output:
% res - [1, 1] - value of Cp criteria
%
% Author: Alexandr Katrutsa, 2016
% E-mail: aleksandr.katrutsa@phystech.edu


p = size(w, 2);
m = size(X, 1);
RSS = par.rss;
res = sumsqr(y - X * w) / RSS - m + 2 * p;
end