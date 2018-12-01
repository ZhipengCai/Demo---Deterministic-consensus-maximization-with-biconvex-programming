% Solve the convex quadratic programming using Frank-Wolfe method
% min_x x^T * Q * x + p^T * x
% subject to A*x + b > 0

function [x] = FrankWolfeQuadProg(Q, p, A, b, tol, x0)
    
    if (nargin < 5); tol = 1e-9; end
    
    
    
end