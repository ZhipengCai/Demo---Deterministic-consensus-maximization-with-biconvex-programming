function [M, inliers, trialcount, outSampleSet] = myHomographyRansac(x, fittingfn, distfn, degenfn, s, t,...
    method, inpSampleSet, feedback, ...
    maxDataTrials, maxTrials,...
    loransacInlierThreshold, matchingScores)

rng shuffle;
if nargin < 7; method='RANSAC'; end;
if nargin < 8; inpSampleSet = []; end;
if nargin < 11; maxTrials = 100000000000;    end;
if nargin < 10; maxDataTrials = 100; end;
if nargin < 9; feedback = 0;        end;
if nargin < 12; loransacInlierThreshold = 0; end;
if nargin < 13; matchingScores = []; end;


[rows, npts] = size(x);

p = 0.99;         % Desired probability of choosing at least one sample
% free from outliers

bestM = NaN;      % Sentinel value allowing detection of solution failure.
trialcount = 0;
bestscore =  0;
N = 1;            % Dummy initialisation for number of trials.
outSampleSet = [];
sampleSetIdx = 0;

%Parameters for PROSAC
n = size(x,2);
psize = s;
sample_num = psize;
M = 10;
Tprime  = 1;
T = M/nchoosek(n,psize);
m = 0;
[~, qinx] = sort(matchingScores);


while (N > trialcount || trialcount <=100)
    
    degenerate = 1;
    count = 1;
    
    % Perform sampling
    while degenerate
        
        if (isempty(inpSampleSet))
            if (isempty(matchingScores)) % If matching score is not provided, just use random sample
                ind = randomsample(npts, s);
            elseif strcmp(method,'PROSAC') % Prosac sampling
                m  = m + 1;
                [pinx, sample_num, T, Tprime] = prosacSample(n, m, qinx, sample_num, s, T, Tprime);
                ind = pinx;
            elseif strcmp(method, 'GMLE') % Guided MLE sampling
                ind = randsample(n,s,true,matchingScores);
            else
                error('Wrong method');
            end
            
        elseif (sampleSetIdx == size(inpSampleSet,1))
            break;
        else
            sampleSetIdx = sampleSetIdx + 1;
            ind = inpSampleSet(sampleSetIdx,:)';
        end
        
        outSampleSet = [outSampleSet; ind']; %save sample to outSampleSet
        
        % Test that these points are not a degenerate configuration.
        degenerate = feval(degenfn, x(:,ind));
        
        if ~degenerate
            
            M = feval(fittingfn, x(:,ind));
            if isempty(M)
                degenerate = 1;
            end
        end
        
        % Safeguard against being stuck in this loop forever
        count = count + 1;
        if count > maxDataTrials
            warning('Unable to select a nondegenerate data set');
            break
        end
    end
    
    % Now, a sample is ready, use it to estimate the model
    [inliers, M] = feval(distfn, M, x, t);
    
    % Find the number of inliers to this model.
    ninliers = length(inliers);
    
    if (ninliers > bestscore & ~isnan(M))    % Largest set of inliers so far...
        bestscore = ninliers;  % Record data for this model
        bestinliers = inliers;
        bestM = M;
%        iter = trialcount;
        
         %Perform Local Update
         if (strcmp(method, 'LORANSAC')||strcmp(method,'F-LORANSAC')) %fixing-LORANSAC
            loInliers = bestinliers;
            loIter = 0;
            
            if(strcmp(method,'LORANSAC'))
                loRansacMaxIter = 10;
                NOSample = max(s,floor(length(loInliers)/2));
            else
                loRansacMaxIter = 50;
                NOSample = min(s*7,length(loInliers));
            end
            
            if (ninliers < s); continue; end
            if (trialcount < 50); loIter = loRansacMaxIter; end
            
            while (loIter < loRansacMaxIter)
                loIter = loIter + 1;
                loind = randsample(loInliers, NOSample);
                loTheta = feval(fittingfn, x(:, loind));
                [loUpdatedInliers, loTheta] = feval(distfn, loTheta, x, t);
                if (numel(loUpdatedInliers) > bestscore)
                    disp(['LO-RANSAC updated from ' num2str(bestscore) ' to ' num2str(numel(loUpdatedInliers)) ' at iter ' num2str(trialcount)]);
                    bestscore = numel(loUpdatedInliers);
                    ninliers = bestscore;
                    bestinliers = loUpdatedInliers;
                    bestM = loTheta;
                end
            end
         end
        
        
        % Update estimate of N, the number of trials to ensure we pick,
        % with probability p, a data set with no outliers.
        fracinliers =  ninliers/npts;
        pNoOutliers = 1 -  fracinliers^s;
        pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
        pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
        N = log(1-p)/log(pNoOutliers);
    end
    
    trialcount = trialcount+1;
    if feedback
        fprintf('trial %d out of %d         \r',trialcount, ceil(N));
    end
    
    % Safeguard against being stuck in this loop forever
    if trialcount > maxTrials
            break
    end
end
fprintf('\n');

if ~isnan(bestM)   % We got a solution
    M = bestM;
    inliers = bestinliers;
else
    M = [];
    inliers = [];
    error('ransac was unable to find a useful solution');
end

end





















% RANSAC - Robustly fits a model to data with the RANSAC algorithm
%
% Usage:
%
% [M, inliers] = ransac(x, fittingfn, distfn, degenfn s, t, feedback, ...
%                       maxDataTrials, maxTrials)
%
% Arguments:
%     x         - Data sets to which we are seeking to fit a model M
%                 It is assumed that x is of size [d x Npts]
%                 where d is the dimensionality of the data and Npts is
%                 the number of data points.
%
%     fittingfn - Handle to a function that fits a model to s
%                 data from x.  It is assumed that the function is of the
%                 form:
%                    M = fittingfn(x)
%                 Note it is possible that the fitting function can return
%                 multiple models (for example up to 3 fundamental matrices
%                 can be fitted to 7 matched points).  In this case it is
%                 assumed that the fitting function returns a cell array of
%                 models.
%                 If this function cannot fit a model it should return M as
%                 an empty matrix.
%
%     distfn    - Handle to a function that evaluates the
%                 distances from the model to data x.
%                 It is assumed that the function is of the form:
%                    [inliers, M] = distfn(M, x, t)
%                 This function must evaluate the distances between points
%                 and the model returning the indices of elements in x that
%                 are inliers, that is, the points that are within distance
%                 't' of the model.  Additionally, if M is a cell array of
%                 possible models 'distfn' will return the model that has the
%                 most inliers.  If there is only one model this function
%                 must still copy the model to the output.  After this call M
%                 will be a non-cell object representing only one model.
%
%     degenfn   - Handle to a function that determines whether a
%                 set of datapoints will produce a degenerate model.
%                 This is used to discard random samples that do not
%                 result in useful models.
%                 It is assumed that degenfn is a boolean function of
%                 the form:
%                    r = degenfn(x)
%                 It may be that you cannot devise a test for degeneracy in
%                 which case you should write a dummy function that always
%                 returns a value of 1 (true) and rely on 'fittingfn' to return
%                 an empty model should the data set be degenerate.
%
%     s         - The minimum number of samples from x required by
%                 fittingfn to fit a model.
%
%     t         - The distance threshold between a data point and the model
%                 used to decide whether the point is an inlier or not.
%
%     feedback  - An optional flag 0/1. If set to one the trial count and the
%                 estimated total number of trials required is printed out at
%                 each step.  Defaults to 0.
%
%     maxDataTrials - Maximum number of attempts to select a non-degenerate
%                     data set. This parameter is optional and defaults to 100.
%
%     maxTrials - Maximum number of iterations. This parameter is optional and
%                 defaults to 1000.
%
% Returns:
%     M         - The model having the greatest number of inliers.
%     inliers   - An array of indices of the elements of x that were
%                 the inliers for the best model.
%
% For an example of the use of this function see RANSACFITHOMOGRAPHY or
% RANSACFITPLANE

% References:
%    M.A. Fishler and  R.C. Boles. "Random sample concensus: A paradigm
%    for model fitting with applications to image analysis and automated
%    cartography". Comm. Assoc. Comp, Mach., Vol 24, No 6, pp 381-395, 1981
%
%    Richard Hartley and Andrew Zisserman. "Multiple View Geometry in
%    Computer Vision". pp 101-113. Cambridge University Press, 2001

% Copyright (c) 2003-2006 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% May      2003 - Original version
% February 2004 - Tidied up.
% August   2005 - Specification of distfn changed to allow model fitter to
%                 return multiple models from which the best must be selected
% Sept     2006 - Random selection of data points changed to ensure duplicate
%                 points are not selected.
% February 2007 - Jordi Ferrer: Arranged warning printout.
%                               Allow maximum trials as optional parameters.
%                               Patch the problem when non-generated data
%                               set is not given in the first iteration.
% August   2008 - 'feedback' parameter restored to argument list and other
%                 breaks in code introduced in last update fixed.
%
