function [ maxVif ] = Vif( X, y, w, par )
% Vif calculates VIF for all features from matrix X and output 
% the maximum value
%   
% [ maxVif ] = calcVif( X )
%
% Input:
% X - [m, n] - variables matrix 
%
% Output:
% maxVif - [1, 1] - maximum Vif along all features
%
% Author: Alexandr Katrutsa, 2016
% E-mail: aleksandr.katrutsa@phystech.edu

if(isempty(X))
    maxVif = Inf;
    return
end

[~, m] = size(X);
VifList = zeros(m,1);

for j = 1:m
    xJ = X(:,j);
    xWithoutJ = X(:,[1:(j-1),(j+1):end]);
    opts.SYM = true;
    regrJ = linsolve((xWithoutJ'*xWithoutJ), (xWithoutJ)'*xJ, opts);
    xHatJ = xWithoutJ*regrJ;% calculated responses for j-th variable
    
    % calculate VIF for j-th variable 
    VifList(j) = sum((xJ-mean(xJ)).^2) / sum((xHatJ - xJ).^2);
end    
maxVif = max(VifList);
end