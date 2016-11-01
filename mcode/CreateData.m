function [ X ] = CreateData( m, features, par )
% Function creates design matrix according to values of the
% features structure fields
%
% Input:
% m - [1, 1] - number objects, number rows in created matrix
% features - structure with fields:
%       rand_features - [1, 1] - number of randon features
%       ortfeat_features - [1, 1] - number of orthogonal features
%       coltarget_features - [1, 1] - number of features collinearing 
%                                   with target vector
%       colfeat_features - [1, 1] - number of features corellating
%                                   with other orthogonal features
%       ortcol_features - [1, 1] - number of features orthogonal 
%                                  target vector and collinearing to each other
% par - structure with fields:
%       multpar - [1, 1] - parameter of multicollinearity, if multpar = 1, 
%                          then all selected features are collinearing
%       target - [m, 1] - target vector
%
% Output:
% X - [m, total_features] - matrix with test data set
%
% Author: Alexandr Katrutsa, 2014-2015 
% E-mail: aleksandr.katrutsa@phystech.edu

random = features.rand_features;
ortfeat = features.ortfeat_features;
coltarget = features.coltarget_features;
colfeat = features.colfeat_features;
ortcol = features.ortcol_features;
k = par.multpar;
target = par.target;

total_features = random + ortfeat + coltarget + colfeat + ortcol;
assert(m >= total_features, 'Not enough objects, objects must be more than features');
assert(ortfeat ~= 1, 'Ortogonal features must be more than 1');
% Create for every kind of features a matrix and concatenate it

%Create random features
mat_random = [];
if random > 0
    mat_random = [rand(m, random - 1), target + 0.01.*randn(m, 1)];
end
% Create orthogonal features, a linear combination of them equals target
% vector
% vec_1 = [target(1:floor(size(target, 1) / 2)); zeros(size(target, 1) - floor(size(target, 1) / 2), 1)];
% vec_2 = [zeros(floor(size(target, 1) / 2), 1); target(1:(size(target, 1) - floor(size(target, 1) / 2)))];

if (ortfeat > 0)
    vec_1 = zeros(size(target, 1), 1);
    vec_1(1:2:end) = target(1:2:end);
    vec_2 = zeros(size(target, 1), 1);
    vec_2(2:2:end) = target(2:2:end);
    if(ortfeat < 3)
        mat_ortfeat = [vec_1, vec_2];
    else
        mat_ort_ortfeat = null([vec_1'; vec_2']);
        mat_ortfeat = [vec_1, vec_2, mat_ort_ortfeat(:, randperm(size(mat_ort_ortfeat, 2), ortfeat-2))];
    end
else
    mat_ortfeat = [];
end

% Create features, which are collinearing to a target vector according to k
mat_coltarget = zeros(m, coltarget);
if coltarget > 0
    mat_ort_coltarget = null(target');
    for i = 1:coltarget
       mat_coltarget(:, i) = k * target + (1 - k) * mat_ort_coltarget(:, i);    
    end
end

% Create features, which are correlated to the orthogonal features 
% from mat_ortfeat 
mat_colfeat = zeros(m, colfeat);
if ortfeat > 1 && colfeat > 0
    idx_first = 1;
    idx_last = 0;
    colfeat_per_ortfeat = floor(colfeat / ortfeat).*ones(1, ortfeat);
    for i = 1:(colfeat - floor(colfeat / ortfeat) * ortfeat)
        colfeat_per_ortfeat(i) = colfeat_per_ortfeat(i) + 1;
    end
    for i = 1:ortfeat
        mat_ort_ortfeat = null(mat_ortfeat(:, i)');
        mat_col_perfeat = zeros(m,colfeat_per_ortfeat(i));
        idx_last = idx_last + colfeat_per_ortfeat(i); 
        for j = 1:colfeat_per_ortfeat(i)
            mat_col_perfeat(:, j) = k * mat_ortfeat(:, j) + (1 - k) * mat_ort_ortfeat(:, j);
        end
        mat_colfeat(:, idx_first:idx_last) = mat_col_perfeat;
        idx_first = idx_last + 1;
    end
end
% Create features, which are orthogonal to the target vector and
% collinearing each other
mat_ortcol = zeros(m, ortcol);
if ortcol > 0
    mat_ort_coltarget = null(target');
    mid = floor(size(mat_ort_coltarget, 2) / 2);
    for i = 1:ortcol
       mat_ortcol(:, i) = k * mat_ort_coltarget(:, mid)  + (1 - k) * mat_ort_coltarget(:, i);   
    end
end
    X = [mat_ortcol, mat_ortfeat, mat_colfeat, mat_coltarget, mat_random];
    %X = X(:, randperm(size(X, 2)));
end
