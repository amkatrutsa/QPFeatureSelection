function [x] = SolveOptProblem(Q, b)
% Function solves the quadratic optimization problem stated to select
% significance and noncollinear features
%
% Input:
% Q - [n, n] - matrix of features similarities
% b - [n, 1] - vector of feature relevances
%
% Output:
% x - [n, 1] - solution of the quadratic optimization problem
%
% Author: Alexandr Katrutsa, 2016 
% E-mail: aleksandr.katrutsa@phystech.edu

[n, ~] = size(Q);
cvx_solver Mosek
cvx_begin
    variable x(n, 1) nonnegative;
    minimize (x'*Q*x - b'*x)
    subject to
        norm(x, 1) <= 1;
cvx_end
end

