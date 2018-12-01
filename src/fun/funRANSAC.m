%% Perform RANSAC for linear fitting
% Method determines the updating method: LO-RANSAC, 
function [theta, inliers, trialcount, sampleSet] = linearFitRANSAC(x, y, th, method, inpSampleSet, loransacInlierThreshold, matchingScores)
    
    if nargin < 4; method = 'RANSAC'; end
    if nargin < 5; inpSampleSet = []; end
    if nargin < 6; loransacInlierThreshold = 0; end;
    if nargin < 7; matchingScores = []; end;
    
    X = [x, y]'; 
    fittingfn = @myFitTchebycheff;
    distfn    = @lineptdist;
    degenfn   = @isdegenerate;

    %maxTrials = 10000;
    maxTrials = 1e3;
    maxDataTrials = 1000;
    feedback = 0;
    
    s = size(x, 2); 
    
    [theta, inliers, trialcount, sampleSet] = myRANSACFun(X, fittingfn, ...
                                                   distfn, degenfn, ...
                                                   s, th, method, ...
                                                   inpSampleSet, feedback, ...
                                                   maxDataTrials,maxTrials,loransacInlierThreshold,...
                                                   matchingScores);
    %compute previous inliers
    theta0 = theta;
    inls0 = find(abs(y-x*theta0)<=th);
    M = ((x(inls0,:)'*x(inls0,:))^-1)*x(inls0,:)'*(y(inls0));
    F = [M(1:3)';M(4:6)';[M(7),1,M(8)]];
    FPrime = enforceFundamentalConstraint(F);
    M = [FPrime(1,:) FPrime(2,:) FPrime(3,1) FPrime(3,3)]';
    theta = M(:);
    inliers = find(abs(y-x*theta)<=th);
    if(numel(inls0)>numel(inliers))
        inliers = inls0;
        theta = theta0;
    end
    
end


function [inliers, xn] = lineptdist(xn, XY, t)

    X = XY(1:end-1, :)';  
    Y = XY(end, :)'; 
    d = X*xn - Y; 
    
    inliers = find(abs(d) < t);
    
end

function r = isdegenerate(X)
    r = 0; 
end

