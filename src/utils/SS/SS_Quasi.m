function [theta, runtime] = SS_Quasi(A, b, c, d, theta0, th, config)
tic
%initialization
if (length(theta0)==9)
    theta0 = theta0./theta0(end);
    theta0 = theta0(1:end-1);
end
theta = theta0;
[~, s0, ~] = compute_residuals_l2(A,b,c,d,theta,th);
s = max(0,s0);
gamma = config.gammaSS;
w = 1./(gamma+s);
%initialize gurobi model for constrains in socp problem (only need to reset model.obj in BCO_l2 function)
NoPoints = numel(d);
d1 = numel(b)/NoPoints;
d2 = size(A,2);
model = genModelSocp_BCO_sedumi(A,b,c,d,th,NoPoints,d1,d2);
maxIter = 1e5;
currF = inf;
for j = 1:maxIter
    preF = currF;
    %solve one LP
    fs = [zeros(length(theta),1); -w]; %obj
    pars.fid=0;
    [~,ys,info] = sedumi(model.At,fs,model.ct,model.K,pars);
    currF = -fs'*ys;
    theta = ys(1:length(theta));

    s = ys(length(theta)+1:end);
    dif = preF-currF;
    w = 1./(gamma+s);
    if(dif<config.QThresh)
        break;
    end
end
runtime = toc;
end
