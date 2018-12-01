% Perform Lo-RANSAC_PR fit for linear data
% x: input data.
% fittingfn: Function to compute a model based on sampled points.
% distfun: Function to compute residuals and inlier w.r.t a specific theta
% s: minimum number of datapoints to fit a model
% t: inlier threshold (epsilon)
% maxDataTrial: Stop RANSAC after maxDataTrials trials. In this experiment, this
%               number was set to a very large number to favor standard stopping
%               criterion with probability of p = 0.99
% maxTrial: Try maximum of maxTrial samples until a non-degenerate sample point set
%           is found.

% M: best model found by Lo-RANSAC

function [M] = myLoRANSAC_PR(x, fittingfn, distfn, degenfn, s, t, maxDataTrials, maxTrials, config)
     
    if nargin < 7; maxDataTrials = 1000; end;
    if nargin < 8; maxTrials = 1e9;    end;
    
    feedback = 0;     % Debug flag. Set this to 1 to print execution steps
    
    [~, npts] = size(x); % Number of points from input data
    
    p = 0.99;         % Desired probability of choosing at least one sample
                      % free from outliers (probably should be a parameter)

    bestM = NaN;      % Sentinel value allowing detection of solution failure.
    trialcount = 0;
    bestscore =  0;
    
    N = 1;            % Dummy initialisation for number of trials.
    
    while (N > trialcount || trialcount < 500) % Force RANSAC to test at least 500 samples before it can quit
        
        % Select at random s datapoints to form a trial model, M.
        % In selecting these points we have to check that they are not in
        % a degenerate configuration.
        degenerate = 1;
        count = 1;
        while degenerate
            ind = randomsample(npts, s);                
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
                                                       
            %bestM = M;             % best model sofar
            
            %run parametric reformulation method
            [bestM, inliers, ninliers] = myParaRef(x(1:(end-1),:)',x(end,:)',t, M, config);
            %compute inliers
            bestscore = ninliers;  % Record data for this model    
            
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
            fprintf('trial %d out of %d    \r',trialcount, ceil(N));
        end

        % Safeguard against being stuck in this loop forever
        if trialcount >= maxTrials
            disp(['RANSAC reaches the maximum number of ' num2str(maxTrials) 'iterations']);
            break;
        end
                
    end    
    if feedback, fprintf('\n'); end
    
    if ~isnan(bestM)   % We got a solution
        M = bestM;       
    else
        M = [];
        error('ransac was unable to find a useful solution');
    end
end
    
