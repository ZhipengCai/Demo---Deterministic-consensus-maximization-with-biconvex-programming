function [theta, inls, iter, runtime, sampleSet] = linearFitFundamental(x,y, th, method, theta0, inpSampleSet, matchingScores,config)

N = size(x,1);
d = size(x,2);
if (nargin < 5);  theta0 = rand(size(x,2),1);  end
if (nargin < 6);  inpSampleSet = []; end;
if (nargin < 7); matchingScores = []; end;

%% Perform Fitting
if  strcmp(method, 'RANSAC') || strcmp(method, 'RANSAC1')
    %% RANSAC
    tic
    [rtheta, inls, iter, sampleSet] = funRANSAC(x, y, th);
    runtime = toc;
elseif strcmp(method, 'LORANSAC')
    %%  LO-RANSAC
    tic
    [rtheta, inls, iter, sampleSet] = funRANSAC(x, y, th, 'LORANSAC', inpSampleSet);
    runtime = toc;
elseif strcmp(method,'FLORANSAC')
    tic
    [rtheta, inls, iter, sampleSet] = funRANSAC(x, y, th, 'F-LORANSAC', inpSampleSet, 0);
    runtime = toc;    
elseif strcmp(method, 'PROSAC') || strcmp(method, 'GMLE')
    tic
    [rtheta, inls, iter, sampleSet] = funRANSAC(x, y, th, method, inpSampleSet, 0, matchingScores); %Set LO-RANSAC Threshold to 100
    runtime = toc;
elseif strcmp(method, 'EP')
    [A, b]=genLinearMaxFSMatrix(x, y, th);
    [x0, y0] = genStartingPoint(A, b, theta0);
    
    alpha = config.alpha;
    kappa = config.kappa;
    
    tic
    while true %(alpha < alphaMax)
        
        [x0, y0, theta, P1, F, Q] = fw(A, b, x0, y0, alpha,config.solver.LP);     
        inls = find(abs(x*theta - y) <=th);
        disp(['Alpha = ' num2str(alpha) ' P=' num2str(P1) ' F=' num2str(F) ' Q=' num2str(Q) ' Inliers = ' num2str(numel(inls))]);
        disp('---------------------------------------------------------------------');
        
        alpha = alpha*kappa;
        if ((Q)<=config.QThresh || alpha > config.maxAlpha) % || F == prevF)
            disp('========================================================');
            break;
        end
    end
    runtime = toc;
    
    rtheta = theta;
    
    iter = -1;
    sampleSet = inpSampleSet;
    
    s = zeros(N,1);
elseif (strcmp(method,'SS'))
    [rtheta, ~, ~, runtime] = SS(x, y, th, theta0, config);
    iter = -1;
    sampleSet = inpSampleSet;
elseif (strcmp(method, 'IBCO'))
    [A, b, c, d] = genQuasiconvexMatrixLinear(x, y);
    [rtheta, runtime] = IBCO_v2(A,b,c,d,theta0,th,config,'fun');
    iter = -1;
    sampleSet = inpSampleSet;
end

theta = rtheta;


%% Find Inliers
if (isempty(theta))
    theta = theta0;
end

F = [theta(1:3)';theta(4:6)';[theta(7),1,theta(8)]];
FPrime = enforceFundamentalConstraint(F);
theta = [FPrime(1,:) FPrime(2,:) FPrime(3,1) FPrime(3,3)]';

inls = find(abs(y-x*theta)<=th);
end











