% linear fitting algorithms for synthetic data
% ---------------inputs-------------------------------------
% 1-2) x, y: data;
% 3) th: inlier threshold
% 4) method: 
%       1. 'RANSAC': ransac
%       2. 'LRS': Lo-RANSAC
%       3. 'FLRS': Fixing Lo-RANSAC
%       4. 'EP': the Exact Penalty method
%       5. 'SS': the Smooth Surrogate method
%       5. 'IBCO': our method
% 5) theta0: initial model
% 6) config: variable for parameter settings
% 7) inpSampleSet: sample set used in RANSAC, used in Lo-RANSAC and
% Fixing Lo-RANSAC
%------------------------------------------------------------
%----------------outputs-------------------------------------
%1) theta: final model
%2) nInls: consensus (number of inliers found)
%3) runtime: algorithm running time
%4) sampleSet: sample set used in RANSAC (set to '[]' for other methods)
%------------------------------------------------------------


function [theta, nInls, runtime, sampleSet] = linearFit(x, y, th, method, theta0, config, inpSampleSet)
if (nargin < 7);  inpSampleSet = []; end;

if strcmp(method, 'RANSAC')   %% RANSAC
    tic
    [theta, ~, ~, sampleSet] = linearFitRANSAC(x, y, th);
    runtime = toc;
elseif strcmp(method, 'LRS')   %% LO-RANSAC
    tic
    [theta] = linearFitRANSAC(x, y, th, 'LORANSAC', inpSampleSet);
    runtime = toc;
elseif strcmp(method, 'FLRS')   %% LO-RANSAC
    tic
    [theta] = linearFitRANSAC(x, y, th, 'F-LORANSAC', inpSampleSet);
    runtime = toc;
elseif (strcmp(method, 'EP'))   %% Exact penalty
    
    [A, b]= genLinearMaxFSMatrix(x, y, th);         % Collect data into matrix A and vector b
    % to form set of constrains Ax <= b
    [x0, y0] = genStartingPoint(A, b, theta0); % From theta0, compute u0, s0, v0
    
    alpha = config.alpha;                      % Inital alpha
    kappa = config.kappa;                      % Incremental rate
    tic
    
    iter = 0;
    while true
        % Execute Frank-Wolfe
        [x0, y0, theta, P1, F, Q] = fw(A, b, x0, y0, alpha, config.solver.LP);          
        if (Q<=config.QThresh || alpha > config.maxAlpha)                               % Q reaches zero, stop
            disp(['EP TERMINATED  as Q(z) reaches ' num2str(Q) ]);
            break;
        end
        % Increase alpha
        alpha = alpha*kappa;
        
    end
    runtime = toc;
elseif (strcmp(method, 'SS'))
    [theta, ~, ~, runtime] = SS(x, y, th, theta0, config);
elseif(strcmp(method, 'IBCO')) 
    [A, b, c, d] = genQuasiconvexMatrixLinear(x, y);
    [theta, runtime] = IBCO_v2(A,b,c,d, theta0,th, config);
end


%% Find Inliers
if (isempty(theta))
    theta = theta0;
end
inls = find(abs(y-x*theta)<=th);
nInls = numel(inls);


end












