function [ w ] = Stepwise(X, y, parameters)
% Function wraps the standard stepwisefit function
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

w = stepwisefit(X, y, 'display', 'off');
end

