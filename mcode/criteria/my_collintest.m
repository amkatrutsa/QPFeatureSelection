function [ VarDecomp ] = my_collintest(X)
% Function computes matrix of the coefficient variance proportion according
% to the Belsley approach to the multicollinear features detection.
% 
% Input:
% X - [m, n] - design matrix
% 
% Output:
% VarDecomp - [n, n] - matrix of the coefficients variance proportion
% 
% Author: Alexandr Katrutsa, 2016
% E-mail: aleksandr.katrutsa@phystech.edu
% 
% Reference: 
% Belsley, D. A., E. Kuh, and R. E. Welsh. Regression Diagnostics. 
% New York, NY: John Wiley & Sons, Inc., 1980. 

[numObs,numVars] = size(X);

% Scale columns to length 1:
colNormsX = sqrt(sum(X.^2));
colNormsX(colNormsX == 0) = 1; % Avoid divide by 0
XS = X./repmat(colNormsX,numObs,1); % Scaled X

% Compute SVD:
[~,S,V] = svd(XS,0);
sValue = diag(S);

% Compute condition indices:
sValue = [sValue; Inf*ones(max(size(S)) - min(size(S)), 1)];
sValue(sValue < eps) = eps; % Avoid divide by 0

% Compute variance decomposition proportions:

PHI = (V.^2)./repmat((sValue.^2)',numVars,1);
phi = sum(PHI,2);
VarDecomp = PHI'./repmat(phi',numVars,1);
end