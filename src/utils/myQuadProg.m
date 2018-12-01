% Solve the convex quadratic programming using Gurobi
% min_x x^T * Q * x + p^T * x
% subject to A*x + b >= 0
%              x >= 0    

function [x] = myQuadProg(Q, p, A, b, tol, x0, sign, method)
if (nargin == 7)
    x = gurobiQuadProg(Q, p, A, b, tol, x0, sign);
    return;
else 
    options = optimoptions('quadprog','Algorithm','interior-point-convex',...
    'Display','off');
    x = quadprog(2.*Q, p, -A, b, [],[], zeros(length(p),1), [], x0, options);
end


end