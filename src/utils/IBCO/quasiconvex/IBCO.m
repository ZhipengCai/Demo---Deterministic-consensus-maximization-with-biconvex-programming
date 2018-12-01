function [theta, runtime] = IBCO(A, b, c, d, theta0, th, config,problem)
if nargin < 8
    problem = 'general';
end
% Normalize matrix if a 3x3 matrix was supplied
% if (length(theta0)==9)
%     theta0 = theta0./theta0(end);
%     theta0 = theta0(1:end-1);
% end

% From theta0, compute u0,s0, v0
%finds the number of inliers
theta = theta0;
NoPoints = numel(d);
d1 = numel(b)/NoPoints;
d2 = size(A,2);

%initialize gurobi model for constrains in socp problem for sx (only need to reset model.obj in BCO_l2 function)
epsilon_ac = th-config.QThresh; %tighten epsilon by accuracy to remove the effect of numerical issue
modelsx = genModelSocp_BCO_sedumi(A,b,c,d,epsilon_ac,NoPoints,d1,d2);

tic

%step 0:
[resn,s,inliers] = compute_residuals_l2(A,b,c,d,theta,th);
u0 = ones(NoPoints,1);
u0(inliers) = 0;
sx0 = [theta; s];

deltaMax = NoPoints-numel(inliers);
deltaMin = 0;

[sx1,u1, theta0, fDeltaMin, fDeltaMinProj] = BCO_l2(A, b, c, d, modelsx, sx0, u0, 0, config, th, epsilon_ac,NoPoints,d1,d2, 'sedumi',problem);


%Find Inliers and exit if problem solved for 0 outlier situation
if fDeltaMinProj <= config.QThresh
    theta = theta0;
    runtime = toc;
    return;
end

nInls = inlierCountQuasi_l2(A,b,c,d,theta0,th);

if deltaMax>(NoPoints-nInls)
    disp(['noInliers: ' num2str(nInls)]);
    deltaMax = NoPoints-nInls;
    theta = theta0;
    sx0=sx1;
    u0=u1;
end

delta = floor(2/3*deltaMax);

while deltaMax>(deltaMin+1)
    update = 0;
    %step 1:
    %execute active-set algorithm
    [sx1, u1, theta0, fDelta, fDeltaProj] = BCO_l2(A, b, c, d, modelsx, sx0, u0, delta, config, th, epsilon_ac,NoPoints,d1,d2, 'sedumi',problem);
    %step 2:
    nInls = inlierCountQuasi_l2(A,b,c,d,theta0,th);
    
    if deltaMax>(NoPoints-nInls)
        disp(['noInliers: ' num2str(nInls)]);
        deltaMax = NoPoints-nInls;
        theta = theta0;
        sx0=sx1;
        u0=u1;
        update = 1;
    end
    %step 3:
    if fDelta <= config.QThresh && update == 1
        delta = floor((deltaMin+deltaMax)/2);
    else
        p = floor(delta - fDelta*(delta-deltaMin)/(fDelta-fDeltaMin));
        deltaMin = delta;
        fDeltaMin = fDelta;
        if (p>deltaMin)&&(p<deltaMax)
            delta = p;
        else
            delta = floor((deltaMin+deltaMax)/2);
        end
    end
end
runtime = toc;
end
