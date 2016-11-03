function [matAlgCrit] = TestAlgCrit(W, crit, X_test, y_test, param)
matAlgCrit = zeros(size(W, 2), length(crit));
beta = lscov(X_test, y_test);
param.rss = sumsqr(y_test - X_test * beta);
for i=1:length(crit)
    for j = 1:size(W, 2)
        X_unnorm = X_test;
        idx_del = W(:, j) == 0;
        X_unnorm(:, idx_del) = [];
        param.X_unnorm = X_unnorm;
        matAlgCrit(j, i) = feval(crit{i}, X_test(:, ~idx_del), y_test, W(~idx_del, j), param);
    end
end
end

