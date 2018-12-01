% Solve the convex quadratic programming using Gurobi
% min_x x^T * Q * x + p^T * x
% subject to A*x + b >= 0
%              x >= 0    
function [x, flag, results] = gurobiQuadProg(Q, p, A, b, tol, x0, sign)
    
    if (nargin < 5); tol = 1e-9; end
    if (nargin < 7); sign = '>'; end
    flag = 0;
        
    clear model;
    model.Q = sparse(Q);
    model.obj = p;
    model.A = sparse(A);
    model.rhs = -b;
    model.sense = sign;
    model.start = x0;
    model.lb = zeros(size(x0));    
    
    params.outflag = 0;
    params.FeasibilityTol = tol;
    params.IntFeasTol = tol;
    params.OptimalityTol = tol;
    params.MarkowitzTol=1e-4;
    params.NormAdjust = 1;
    params.ScaleFlag = 1;
    params.PerturbValue = 0;
%     params.MIPGapAbs = 0;
  %   params.Sifting = 2;
    %params.SiftMethod = 1;
%     params.MIPGap = 0;
    params.Quad = 1;
    
    params.OutputFlag = 0;
    params.Method  = -1;
    %params.BarIterLimit = 50000;
    
    results = gurobi(model, params);
    x = [];
    if (~(strcmp(results.status,'OPTIMAL') || strcmp(results.status,'SUBOPTIMAL')))
        flag = -1;
    else 
        x = results.x;
    end
    
%     
%     clear model;
% names = {'x', 'y', 'z'};
% model.varnames = names;
% model.Q = sparse([1 0.5 0; 0.5 1 0.5; 0 0.5 1]);
% model.A = sparse([1 2 3; 1 1 0]);
% model.obj = [2 0 0];
% model.rhs = [4 1];
% model.sense = '>';
% 
% gurobi_write(model, 'qp.lp');
% 
% results = gurobi(model);
% 
% for v=1:length(names)
%     fprintf('%s %e\n', names{v}, results.x(v));
% end
% 
% fprintf('Obj: %e\n', results.objval);
% 
% model.vtype = 'B';
% 
% results  = gurobi(model);
% 
% for v=1:length(names)
%     fprintf('%s %e\n', names{v}, results.x(v));
% end
% 
% fprintf('Obj: %e\n', results.objval);

    
end