clear all
addpath('./criteria')
addpath('./alg')
addpath('../data/');
addpath('./mi/');
addpath('./mi2/');
%% Parameters for data sets generation 
rng(0);
% Number of samples
objects = 1000;
% Number of random features
features.rand_features = 50;
% Number of orthogonal features
features.ortfeat_features = 0;
% Number of features collinearing with target vector
features.coltarget_features = 0;
% Number of features corellating with other orthogonal features
features.colfeat_features = 0;
% Number of features orthogonal to target vector and collinearing to each
% other ones
features.ortcol_features = 0;
% A parameter of multicollinearity
param.multpar = 0.8;
% A set of the considered feature selection methods
alg = {'Lasso', 'LARS', 'Stepwise', 'ElasticNet', 'Ridge', 'Genetic'};
% alg = {'Lasso', 'Ridge', 'ElasticNet'};
% A set of the considered criteria 
% crit = {'stability', 'Cp', 'RSS', 'CondNumber', 'Vif', ...
%          'Rsq_adj', 'bic', 'complexity'};
crit = {'complexity', 'Cp', 'RSS', 'CondNumber', 'Vif', 'bic'};
% Number of the iteration in AlgCrit function
param.iter = 1;
param.crit = crit;
% A limit error
param.s_0 = 0.5;
param.threshold = 10^(-10); % to shrink the small coefficients in w* 
param.data = 'real'; % or 'real' 'artificial'
% Parameters of the real data set
param.real_data_filename = 'BP50GATEST.mat';
param.real_data_X = 'bp50_s1d_ll_a';
param.real_data_y = 'bp50_y1_ll_a';
% Parameters for genetic algorithm
param.Genetic.nGenerations = 10;
param.Genetic.nIndividuals = 20;
% Generate the target vector
param.target = randi(1.5 * objects, objects, 1);
% [X, y] = CreateData2(objects, features, param);
% load(param.real_data_filename);
% X = eval(param.real_data_X);
% y = eval(param.real_data_y);
X = CreateData(objects, features, param);
y = param.target;
% Optional normalization
len = sum(X.^2).^0.5;
X = X./repmat(len, size(X, 1), 1);
y = y ./ norm(y);
%% Slit test and train set
test_set_ratio = 0.7;
X_train = X(1:floor(test_set_ratio*size(X, 1)), :);
y_train = y(1:floor(test_set_ratio*size(X, 1)));
X_test = X(floor(test_set_ratio*size(X, 1)) + 1:size(X, 1), :);
y_test = y(floor(test_set_ratio*size(X, 1)) + 1:size(X, 1));
%% Create and solve quadratic optimization problem
% String indications for way to compute similarities and relevances
sim = 'correl';
rel = 'correl';
[Q, b] = CreateOptProblem(X_train, y_train, sim, rel);
x = SolveOptProblem(Q, b);
%% Tune significance threshold
threshold = sort(x)';
% threshold = 1e-3:1e-5:0.1;
rss = zeros(1, length(threshold));
rss_test = zeros(1, length(threshold));
stability = zeros(1, length(threshold));
vif = zeros(1, length(threshold));
complexity = zeros(1, length(threshold));
BIC = zeros(1, length(threshold));
cp = zeros(1, length(threshold));
A = zeros(length(threshold), size(X, 2));
for i=1:length(threshold)
    fprintf('i = %d\n', i);
    active_idx = x >= threshold(i);
    A(i, :) = active_idx;
    if sum(active_idx) > 0
        lm = fitlm(X_train(:, active_idx), y_train);
        w = lm.Coefficients.Estimate(2:end);
%         w = lscov(X_train(:, active_idx), y_train);
        rss(i) = sumsqr(X_train(:, active_idx)*w - y_train);
        rss_test(i) = sumsqr(X_test(:, active_idx)*w - y_test);
%         stability(i) = CondNumber(X_train(:, active_idx));
        vif(i) = Vif(X_train(:, active_idx));
        complexity(i) = sum(active_idx);
%         BIC(i) = bic(X_train(:, active_idx), y_train, w);
%         p.rss = rss(i);
%         cp(i) = Cp(X_train(:, active_idx), y_train, w, p);
    else
        break;
    end
end
% threshold = threshold(rss > 0);
% rss = rss(rss > 0);
% rss_test = rss_test(rss_test > 0);
% cp = cp(rss > 0);
% vif = vif(rss > 0);
% complexity = complexity(rss > 0);
% stability = stability(rss > 0);
% BIC = BIC(rss > 0);
%% Perform QP feature selection and build model
% [~, idx_min_rss_test] = min(rss_test);
% num_selected_feat = complexity(idx_min_rss_test);
% [sort_x, idx] = sort(x);
% idx = idx(end-num_selected_feat:end);
% % threshold = 1e-3;
% % idx = x > threshold;
% w_qp = lscov(X_train(:, idx), y_train);
% w = zeros(size(X_train, 2), 1);
% w(idx) = w_qp;
%% Function to get tables with values of the evaluation criteria fron all array 'crit'
% [matAlgCrit, W] = AlgCrit(alg, crit, X_train, y_train, param);
% % W(:, end+1) = w;
% crit = {'RSS', 'CondNumber'};
% matAlgCrit_test = TestAlgCrit(W, crit, X_test, y_test, param);