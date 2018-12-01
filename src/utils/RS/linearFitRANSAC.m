%% Perform RANSAC for linear fitting
% Method determines the updating method: LO-RANSAC, 
function [theta, inliers, trialcount, sampleSet] = linearFitRANSAC(x, y, th, method, inpSampleSet, matchingScores)
    
    if nargin < 4; method = 'RANSAC'; end
    if nargin < 5; inpSampleSet = []; end
    if nargin < 6; matchingScores = []; end;
    
    X = [x, y]'; 
    fittingfn = @myFitTchebycheff;
    distfn    = @lineptdist;
    degenfn   = @isdegenerate;

    maxTrials = 11000000;
    maxDataTrials = 1000;
    feedback = 0;
    
    s = size(x, 2); 
    [theta, inliers, trialcount, sampleSet] = myRANSAC(X, fittingfn, ...
                                                   distfn, degenfn, ...
                                                   s, th, method, ...
                                                   inpSampleSet, feedback, ...
                                                   maxDataTrials,maxTrials,...
                                                   matchingScores);
    %perform LS on the final consensus set
    %compute initial inliers
    theta0 = theta;
    inls0 = find(abs(y-x*theta)<=th);
    theta = ((x(inls0,:)'*x(inls0,:))^-1)*x(inls0,:)'*(y(inls0));
    inliers = find(abs(y-x*theta)<=th);
    if (numel(inls0)>inliers)
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

