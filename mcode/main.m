clear all
addpath('./criteria')
addpath('./alg')
addpath('../data/');
addpath('./mi/');
%% Parameters for data sets generation 
% rng(0);
% Number of samples
objects = 1000;
% Number of random features
features.rand_features = 0;
% Number of orthogonal features
features.ortfeat_features = 10;
% Number of features collinearing with target vector
features.coltarget_features = 0;
% Number of features corellating with other orthogonal features
features.colfeat_features = 40;
% Number of features orthogonal to target vector and collinearing to each
% other ones
features.ortcol_features = 0;
% A parameter of multicollinearity
param.multpar = 0.8;
% Generate the target vector
param.target = randi(1.5 * objects, objects, 1);
% A set of the considered feature selection methods
alg = {'Lasso', 'LARS', 'Stepwise', 'ElasticNet', 'Ridge', 'Genetic'};
% A set of the considered criteria 
crit = {'stability', 'Cp', 'RSS', 'CondNumber', 'Vif', ...
         'Rsq_adj', 'bic', 'complexity'};
% Number of the iteration in AlgCrit function
param.iter = 1;
param.crit = crit;
% A limit error
param.s_0 = 0.5;
param.threshold = 10^(-6); % to shrink the small coefficients in w* 
param.data = 'artificial'; % or 'real'
% Parameters of the real data set
param.real_data_filename = 'BP50GATEST.mat';
param.real_data_X = 'bp50_s1d_ll_a';
param.real_data_y = 'bp50_y1_ll_a';
% Parameters for genetic algorithm
param.Genetic.nGenerations = 10;
param.Genetic.nIndividuals = 20;
% parameters.target = randi(1.5 * objects, objects, 1);
parameters.target = randn(objects, 1) * 10;
[X, y] = CreateData(objects, features, param);
% Optional normalization
% len = sum(X.^2).^0.5;
% X = X./repmat(len, size(X, 1), 1);
% y = y ./ norm(y);
%% Create and solve quadratic optimization problem
% String indications for way to compute similarities and relevances
sim = 'correl';
rel = 'signif';
[Q, b] = CreateOptProblem(X, y, sim, rel);
x = SolveOptProblem(Q, b);
%% Tune significance threshold
threshold = sort(x)';
rss = zeros(1, length(threshold));
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
        lm = fitlm(X(:, active_idx), y);
        w = lm.Coefficients.Estimate(2:end);
        rss(i) = norm(X(:, active_idx)*w - y)^2;
        stability(i) = CondNumber(X(:, active_idx));
        vif(i) = Vif(X(:, active_idx));
        complexity(i) = sum(active_idx);
        BIC(i) = bic(X(:, active_idx), y, w);
        p.rss = rss(i);
        cp(i) = Cp(X(:, active_idx), y, w, p);
    else
        break;
    end
end
threshold = threshold(rss > 0);
rss = rss(rss > 0);
cp = cp(rss > 0);
vif = vif(rss > 0);
complexity = complexity(rss > 0);
stability = stability(rss > 0);
BIC = BIC(rss > 0);
%% Function to get tables with values of the evaluation criteria fron all array 'crit'
% matAlgCrit = AlgCrit(alg, crit, objects, features, param);