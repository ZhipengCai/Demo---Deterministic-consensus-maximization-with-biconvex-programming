%initialize gurobi model for constrains in socp problem for sx (only need to reset model.obj in BCO_l2 function)
function modelsx = genModelSocp_BCO(A,b,c,d,epsilon,n,d1,d2)
%sense of the obj
modelsx.modelsense = 'min';
%linear constraints
A1 = [-c', -eye(n), zeros(n,d1*n), eye(n)];
A2 = [-A/epsilon, zeros(n*d1,n), eye(n*d1), zeros(n*d1,n)];
%% first try
% modelsx.A = sparse([A1; A2]);
% modelsx.rhs = [d'; reshape(b/epsilon,n*d1,1)];
% modelsx.sense = '=';
% %quadratic constraints
% QcZeros = zeros(d2+n*(d1+2));
% qZeros = zeros(d2+n*(d1+2),1);
% for i = 1:n
%     Qc = QcZeros;
%     Qc(d2+n+(i-1)*d1+1:d2+n+i*d1,d2+n+(i-1)*d1+1:d2+n+i*d1) = eye(d1);
%     Qc(d2+n+d1*n+i,d2+n+d1*n+i) = -1;
%     modelsx.quadcon(2*i-1).Qc = sparse(Qc);
%     modelsx.quadcon(2*i-1).q = qZeros;
%     modelsx.quadcon(2*i-1).rhs = 0.0;
%     %si>=0
%     modelsx.quadcon(2*i).Qc = sparse(QcZeros);
%     q = qZeros;
%     q(d2+i) = -1;
%     modelsx.quadcon(2*i).q = q;
%     modelsx.quadcon(2*i).rhs = 0.0;
% end
% modelsx.lb = -inf(d2+n+d1*n+n,1);
% modelsx.ub = inf(d2+n+d1*n+n,1);
%% another try
modelsx.A = sparse([A1; A2]);
modelsx.rhs = [d'; reshape(b/epsilon,n*d1,1)];
modelsx.sense = '=';
for i = 1:n
    modelsx.cones(i).index = [d2+n+d1*n+i, (d2+n+d1*(i-1)+1):(d2+n+d1*i)];
end
modelsx.lb = [-inf(d2,1); zeros(n,1); -inf(d1*n,1); zeros(n,1)];
end