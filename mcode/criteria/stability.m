function [ d ] = stability( X, y, w, par )
% Function computes the number of features, 
% after their deleting the error function is less than par.s_0. 
% Deleting features is implemented through the Belsley procedure.     
%   
% Input:
% X - [m, p] - design matrix with shrinkage number of predictors
% y - [m, 1] - target vector
% w - [p, 1] - vector of parameters, getting from algorithm, 
%              which is tested
% par - structure - structure with additional parameters:
%       par.s_0 - [1, 1] - limit acccepted error rate
%       par.X_unnorm - [m, p] - design matrix with shrinkage number of predictors
%                               but not normalized, because of the 
%                               Belsley diagnostic implementation
%
% Output:
% d - [1,1] - maximum number of possibly deleting features
%
% Author: Alexandr Katrutsa, 2016
% E-mail: aleksandr.katrutsa@phystech.edu

if(isempty(X))
    d = Inf;
    return
end
s_0 = par.s_0;
X_unnorm = par.X_unnorm;
idx_all_features = 1:max(size(w));
d = 0;
S = sumsqr(y - X * w);
while (S < s_0) && (size(idx_all_features, 2) > 1)
    idx_del = algBelsley(X_unnorm, idx_all_features);
    idx_all_features(idx_all_features == idx_del) = [];
    d = d + 1;    
    S = sumsqr(y - X(:, idx_all_features) * w(idx_all_features));
end

end

function [ intIdxDelFeature ] = algBelsley( X, idxFeatures )
% Find index the worst feature, the most collinear, through the Belsley
% diagnosis
% [ intIdxDelFeature ] = algBelsley( X, idxFeatures )
% 
% Input:
% X - [m, p] - full design matrix
% idxFeatures - [1, k] - vector containing indices of currently used features
% 
% Output:
% intIdxDelFeature - [1, 1] - index of the deleted feature from idxFeature
%                             vector

VarDecomp = my_collintest(X(:, idxFeatures));
[~, idxMaxVarProp] = max(VarDecomp(end, :));
intIdxDelFeature = idxFeatures(idxMaxVarProp);
end
