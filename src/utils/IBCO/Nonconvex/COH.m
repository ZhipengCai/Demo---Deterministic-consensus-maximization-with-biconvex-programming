% A,b,c,d are coefficients in quasiconvex residuals: ||A_i^T*theta+b_i||/(c_i^T*x+d_i)
%x: d1*1 vector
%A = [A_1, ..., A_n] (d1*(d2*n) matrix)
%b = [b_1, ..., b_n] (d2*n matrix)
%c = [c_1, ..., c_n] (d1*n matrix)
%d = [d_1, ..., d_n] (1*n matrix)
%theta is a d1*1 vector
%constraintType tells the type of convex constraint for theta, if not
%specified, the problem is the same as the one handled by ibco
function [theta, inls, runtime] = COH(theta0, th, config, A, b, constraintType, c, d)
if nargin<6
    constraintType = 'no'; %no constriant
end
if nargin<8
    c = [];
    d = [];
end
%size of the problem
d1 = numel(theta0);
[d2,n] = size(b);
%start algorithm
tic;
%init
theta = theta0;
[resn,s,inls] = compute_residuals_l2(A',b,c,d,theta,th);
%init variable u and sx
u0 = ones(n,1);
u0(inls) = 0;
sx0 = [theta; s];
%init the upper and lower bound of delta
disp(['initial nInls: ', num2str(numel(inls))]);

deltaMax = n-numel(inls);
deltaMin = 0;

disp(['deltaMax = ', num2str(deltaMax), '; deltaMin = ', num2str(deltaMin)]);

%generate model for Yamip
th_accu = th-config.QThresh;%removing the numerical accuracy gap
modelYamip = genModelSDP_COH_Yamip(A,b,c,d,th_accu,n,d1,d2,constraintType);
%solve for delta = 0
[sx1, u1, theta1, fDeltaMin, fDeltaMinProj] = COH_OneStep(A,b,c,d,modelYamip,sx0,u0,0,config,th,th_accu,n,d1,d2,constraintType);

%Find Inliers and exit if problem solved for 0 outlier situation
if fDeltaMinProj <= config.QThresh
    theta = theta1;
    inls = [1:n]';
    runtime = toc;
    return;
end

%calculate inlier set for new solution
[~,~,inls1] = compute_residuals_l2(A',b,c,d,theta1,th);
nInls1 = numel(inls1);
%update solution if better solutions are found
if deltaMax>(n-nInls1)
    disp(['noInliers: ' num2str(nInls1)]);
    deltaMax = n-nInls1;
    theta = theta1;
    inls = inls1;
    sx0=sx1;
    u0=u1;
end
%init delta before iteration
delta = floor(2/3*deltaMax);
disp(['-------------------------------------------------']);

%bisecting the range between deltaMax and deltaMin
while deltaMax>(deltaMin+1)

    disp(['deltaMax = ', num2str(deltaMax), '; deltaMin = ', num2str(deltaMin)]);
    %execute one step of COH algorithm
    [sx1, u1, theta1, fDelta, fDeltaProj] = COH_OneStep(A,b,c,d,modelYamip,sx0,u0,delta,config,th,th_accu,n,d1,d2,constraintType);

    %calculate inlier set for new solution
    [~,~,inls1] = compute_residuals_l2(A',b,c,d,theta1,th);
    nInls1 = numel(inls1);
    
    %update solution if better solutions are found
    if deltaMax>(n-nInls1)
        disp(['noInliers: ' num2str(nInls1)]);
        deltaMax = n-nInls1;
        theta = theta1;
        inls = inls1;
        sx0=sx1;
        u0=u1;
    else
        deltaMin = delta;
    end
    %just do bisection for now, check whether secant method is more
    %effective
    delta = floor((deltaMin+deltaMax)/2);
    disp(['-------------------------------------------------']);
end


runtime = toc;
end