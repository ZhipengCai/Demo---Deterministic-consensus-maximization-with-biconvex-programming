function [M, inliers, iter, outSampleSet] = myRANSAC(x, fittingfn, distfn, degenfn, s, t, method, inpSampleSet, feedback, maxDataTrials, maxTrials, loransacInlierThreshold, matchingScores)

% Test number of parameters
%error ( nargchk ( 6, 9, nargin ) );

if nargin < 7; method='RANSAC'; end;
if nargin < 8; inpSampleSet = []; end;
if nargin < 11; maxTrials = 1000;    end;
if nargin < 10; maxDataTrials = 100; end;
if nargin < 9; feedback = 0;        end;
if nargin < 12; loransacInlierThreshold = 0; end;

[rows, npts] = size(x);


p = 0.99;         % Desired probability of choosing at least one sample
% free from outliers (probably should be a parameter)

bestM = NaN;      % Sentinel value allowing detection of solution failure.
trialcount = 0;
bestscore =  0;
ransacBestScore = 0;
N = 1;            % Dummy initialisation for number of trials.

n = size(x,2);
dim = size(x,1)-1;
psize = s;
sample_num = psize;
M = 1000;
Tprime = 1;
T = M/nchoosek(n, psize);
m = 0;
[sortedScores, qinx] = sort(-matchingScores);

outSampleSet = [];
sampleSetIdx = 0;
while (N > trialcount || trialcount < 100)%||(trialcount <= size(inpSampleSet,1))
    
    
    if (~isempty(inpSampleSet) && (sampleSetIdx == size(inpSampleSet,1)))
        break;
    end
    
    % Select at random s datapoints to form a trial model, M.
    % In selecting these points we have to check that they are not in
    % a degenerate configuration.
    degenerate = 1;
    count = 1;
    while degenerate
        % Generate s random indicies in the range 1..npts
        if (isempty(inpSampleSet))
            
            if (isempty(matchingScores))
                ind = randomsample(npts, s);
            elseif strcmp(method,'PROSAC')
                m  = m + 1;
                [pinx, sample_num, T, Tprime] = prosacSample(n, m, qinx, sample_num, s, T, Tprime);
                ind = pinx;
            elseif strcmp(method, 'GMLE')
                ind = randsample(n,s,true,matchingScores(qinx));
            else
                error('Wrong method');
            end
            
        elseif (sampleSetIdx == size(inpSampleSet,1))
            break;
        else
            sampleSetIdx = sampleSetIdx + 1;
            ind = inpSampleSet(sampleSetIdx,:)';
        end
        
        outSampleSet = [outSampleSet; ind'];
        
        % Test that these points are not a degenerate configuration.
        degenerate = feval(degenfn, x(:,ind));
        
        if ~degenerate
            % Fit model to this random selection of data points.
            % Note that M may represent a set of models that fit the data in
            % this case M will be a cell array of models
            M = feval(fittingfn, x(:,ind));
            
            
            % Depending on your problem it might be that the only way you
            % can determine whether a data set is degenerate or not is to
            % try to fit a model and see if it succeeds.  If it fails we
            % reset degenerate to true.
            
            if isempty(M)||((sum(isinf(M))>=1) || (sum(isnan(M))>=1))
                degenerate = 1;
            end
        end
        
        % Safeguard against being stuck in this loop forever
        count = count + 1;
        if count > maxDataTrials
            warning('Unable to select a nondegenerate data set');
            break;
        end
    end
    
    
    if ((sum(isinf(M))>=1) || (sum(isnan(M))>=1))
        continue;
    end
    
    %enforce the constraint on M for fundamental matrix estimation
    F = [M(1:3)';M(4:6)';[M(7),1,M(8)]];
    FPrime = enforceFundamentalConstraint(F);
    M = [FPrime(1,:) FPrime(2,:) FPrime(3,1) FPrime(3,3)]';
    
    
    % Once we are out here we should have some kind of model...
    % Evaluate distances between points and model returning the indices
    % of elements in x that are inliers.  Additionally, if M is a cell
    % array of possible models 'distfn' will return the model that has
    % the most inliers.  After this call M will be a non-cell object
    % representing only one model.
    [inliers, M] = feval(distfn, M, x, t);
    
    
    % Find the number of inliers to this model.
    ninliers = length(inliers);
    
    if ninliers > bestscore    % Largest set of inliers so far...
        bestscore = ninliers;  % Record data for this model
        ransacBestScore = ninliers;
        
        bestinliers = inliers;
        bestM = M;
        iter=trialcount;
        xdata = x(1:end-1,:)';
        ydata = x(end,:)';
        
        if (strcmp(method,'LORANSAC'))
            loRansacMaxIter = 10;
            loInliers = bestinliers;
            %xInliers = x(bestinliers,:)
            loIter = 0;
            irlsSteps = 10;
            if (ninliers<s); continue;   end
            %if (ninliers < loransacInlierThreshold); loIter = loRansacMaxIter; end;
            
            while (loIter < loRansacMaxIter)
                loIter = loIter + 1;
                % Start sampling from the InlierSet (higher than
                % minimal subsets)
                loind = randsample(loInliers, max(s, length(loInliers)/2));
                % Re-estimate the model
                loTheta = feval(fittingfn, x(:,loind));
                
                %Peform iteratively reweighted least square
                th_multiplier = 3; th_step_size = (th_multiplier*t - t)./irlsSteps;
                for loirls = 1:irlsSteps
                    [loInls, loTheta] = feval(distfn, loTheta, x, t*th_multiplier - loirls*th_step_size);
                    loX = x(1:dim,loInls)'; loY = x(end, loInls)';
                    loRes = abs(loX*loTheta - loY);
                    loWeight = 1./loRes; W = diag(loWeight);
                    %Weighted LS:
                    loTheta = (loX'*W*loX)\(loX'*W*loY);
                end
                
                if ((sum(isinf(loTheta))>=1) || (sum(isnan(loTheta))>=1))
                    continue;
                end
                
                F = [loTheta(1:3)';loTheta(4:6)';[loTheta(7),1,loTheta(8)]];
                FPrime = enforceFundamentalConstraint(F);
                loTheta = [FPrime(1,:) FPrime(2,:) FPrime(3,1) FPrime(3,3)]';
                
                [loUpdatedInliers, loTheta] = feval(distfn, loTheta, x, t);
                if (numel(loUpdatedInliers)>bestscore)
                    disp(['LO-RANSAC updated from ' num2str(bestscore) ' to ' num2str(numel(loUpdatedInliers)) ' at iter ' num2str(trialcount)]);
                    bestscore = numel(loUpdatedInliers);
                    ninliers = bestscore;
                    bestinliers = loUpdatedInliers;
                    bestM = loTheta;
                    iter=trialcount;
                end
                
            end
        elseif(strcmp(method,'F-LORANSAC'))
            loRansacMaxIter = 50;
            loInliers = bestinliers;
            %xInliers = x(bestinliers,:)
            loIter = 0;
            irlsSteps = 10;
            if (ninliers<s); continue;   end
            if (iter < 50); loIter = loRansacMaxIter; end;
            
            while (loIter < loRansacMaxIter)
                loIter = loIter + 1;
                % Start sampling from the InlierSet (higher than
                % minimal subsets)
                loind = randsample(loInliers, min(s*7, length(loInliers)));
                % Re-estimate the model
                loTheta = feval(fittingfn, x(:,loind));
                
                %Peform iteratively reweighted least square
                th_multiplier = 3; th_step_size = (th_multiplier*t - t)./irlsSteps;
                for loirls = 1:irlsSteps
                    [loInls, loTheta] = feval(distfn, loTheta, x, t*th_multiplier - loirls*th_step_size);
                    loX = x(1:dim,loInls)'; loY = x(end, loInls)';
                    loRes = abs(loX*loTheta - loY);
                    loWeight = 1./loRes; W = diag(loWeight);
                    %Weighted LS:
                    loTheta = (loX'*W*loX)\(loX'*W*loY);
                end
                
                
                if ((sum(isinf(loTheta))>=1) || (sum(isnan(loTheta))>=1))
                    continue;
                end
                
                F = [loTheta(1:3)';loTheta(4:6)';[loTheta(7),1,loTheta(8)]];
                FPrime = enforceFundamentalConstraint(F);
                loTheta = [FPrime(1,:) FPrime(2,:) FPrime(3,1) FPrime(3,3)]';
                
                [loUpdatedInliers, loTheta] = feval(distfn, loTheta, x, t);
                if (numel(loUpdatedInliers)>bestscore)
                    disp(['Fixing LO-RANSAC updated from ' num2str(bestscore) ' to ' num2str(numel(loUpdatedInliers)) ' at iter ' num2str(trialcount)]);
                    bestscore = numel(loUpdatedInliers);
                    ninliers = bestscore;
                    bestinliers = loUpdatedInliers;
                    bestM = loTheta;
                    iter=trialcount;
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
    
    
    
    if (strcmp(method, 'LORANSAC') || strcmp(method, 'F-LORANSAC')) % prevent early exit for LORANSAC and F-LORANSAC
        N = trialcount  + 1;
    end
    
    
    if feedback
        fprintf('trial %d out of %d         \r',trialcount, ceil(N));
    end
    
    % Safeguard against being stuck in this loop forever
    if trialcount >= maxTrials
        break;
    end
    
    if (~isempty(inpSampleSet) && (trialcount > size(inpSampleSet,1)))
        break;
    end
end
iter = trialcount;
if feedback, fprintf('\n'); end

if ~isnan(bestM)   % We got a solution
    M = bestM;
    inliers = bestinliers;
else
    M = [];
    inliers = [];
    error('ransac was unable to find a useful solution');
end

