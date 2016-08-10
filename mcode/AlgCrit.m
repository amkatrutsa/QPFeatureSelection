function [ matAlgCrit ] = AlgCrit(alg, crit, objects, features, parameters)
% Function computes values of every considered criteria for result
% given every feature selection method
%
% Input:
% alg - cell array with names of the algorithms to solve regression problem 
%                using feature selection for given design matrix X and target vector y
% crit - cell array with name of criterias using to test feature selection methods from alg
% objects - [1, 1] - number of the objects in generated data set, 
%                    the number of rows in the matrix X
% features - [1, 1] - number of features in generated data set,
%                     the number of columns in the matrix X,
%                     the dimension of the data
% parameters - structure with following parameters:
%              multpar - [1, 1] - parameter of the multicollinearity, 
%                                 if it equals 1, the data is full correlated
%              target - [objects, 1] - the target vector   
%              iter - [1, 1] - number of iteration for average matrix matAlgCrit
%              threshold - [1, 1] - value of the element of w, which equal 0 
%              data - string - indicate what type of data ust be used for
%              testing: 'artificial' or 'real'
%              real_data_filename - string - name of the mat-file with real
%              data set
%              real_data_X - string - name of the variable from the
%              mat-file, which corresponds to the real data design matrix
%              real_data_y - string - name of the variable from the
%              mat-file, which corresponds to the real data target vector
% 
% Output:
% matAlgCrit - [length(alg), length(crit)] - matrix contain the values, 
%                                            which are returned every criteria for 
%                                            every feature selection method 
%
% Author: Alexandr Katrutsa, 2016
% E-mail: aleksandr.katrutsa@phystech.edu

iter = parameters.iter;
threshold = parameters.threshold;
if strcmp(parameters.data, 'real')
    crit(6) = [];
end
matAlgCrit = zeros(length(alg), length(crit));
for it = 1:iter
    if (strcmp('real', parameters.data))
        load(parameters.real_data_filename);
        X = eval(parameters.real_data_X);
        y = eval(parameters.real_data_y);
    elseif (strcmp('artificial', parameters.data))
        parameters.target = randi(1.5 * objects, objects, 1);
        X = CreateData(objects, features, parameters);
        y = parameters.target;
    else
        fprintf('The field parameters.data must be equal "artificial" or "real"!\n')    
    end
    len = sum(X.^2).^0.5;
    X_def = X;
    X = X./repmat(len, size(X, 1), 1);
    y = y ./ norm(y);
    beta = lscov(X, y);
    rss_all = sumsqr(y - X * beta);
    tss_all = sumsqr(y - mean(y));
    W = zeros(size(X, 2), length(alg));
    par.s_0 = parameters.s_0;
    par.rss = rss_all;
    par.tss = tss_all;
    parameters.crit = crit;
    parameters.rss = rss_all;
    for i = 1:length(alg)
        fprintf('Current algorithm is %s...\n', alg{i});
        X_sh = X;
        X_unnorm = X_def;
        parameters.X_unnorm = X_unnorm;
        w = feval(alg{i}, X, y, parameters);
        W(:, i) = w;
        idx_del = abs(w) < threshold;
        w(idx_del) = [];
        W(idx_del, i) = 0;
        X_sh(:, idx_del) = [];
        X_unnorm(:, idx_del) = [];
        par.X_unnorm = X_unnorm;
        for j = 1:length(crit)
            matAlgCrit(i, j) = matAlgCrit(i, j) + feval(crit{j}, X_sh, y, w, par);
        end
    end
end
matAlgCrit = matAlgCrit / iter;
end
