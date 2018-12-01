function [theta, inls, runtime, sampleSet] = imageHomographyFit(x1, x2, th, method, theta0, inpSampleSet, matchingScores,config)


if (nargin<7); matchingScores = []; end

normalizeData = 1;

[A, b, c, d] = genMatrixHomography(x1, x2);

sampleSet = inpSampleSet;

if (strcmp(method,'RANSAC') || strcmp(method, 'LORANSAC') || strcmp(method, 'GREEDY_RANSAC'))
    tic
    [ransacH, rinls, ransacIters, T1, T2, sampleSet] =  ransacfithomography(x1, x2, th,normalizeData, method, inpSampleSet);
    runtime = toc;
    theta = ransacH(:);
elseif (strcmp(method, 'FLORANSAC'))
    tic
    [ransacH, rinls, ransacIters, T1, T2, sampleSet] = ransacfithomography(x1, x2, th,normalizeData, 'F-LORANSAC', inpSampleSet, 0);
    runtime = toc;
    theta = ransacH(:);
elseif (strcmp(method, 'PROSAC') || strcmp(method,'GMLE')) %Provide RASAC with matching scores.
    tic
    [ransacH, rinls, ransacIters, T1, T2, sampleSet] =  ransacfithomography(x1, x2, th,normalizeData, method, inpSampleSet, 0,matchingScores);
    runtime = toc;
    theta = ransacH(:);
elseif (strcmp(method,'EP'))
    [lA, lb] = genLinearMatrixFromQuasiconvex(A, b, c, d, th);
    alpha = config.alpha;
    kappa = config.kappa;
    QThresh = config.QThresh;
    alphaMax = 1e9;
    %Normalize theta0 to make H33 = 1;
    if (length(theta0)==9)
        theta0 = theta0./theta0(end);
        theta0 = theta0(1:end-1);
    end
    
    [x0, y0] = genStartingPoint(lA, lb, theta0);
    
    tic
    while (true)
        [x0, y0, theta, P, F, Q, iters] = MaxFS_LPEC(lA, lb, c, d, th, x0, y0, alpha, 0,config.solver.LP);
        
        %theta(end) = 1;
        
        [res1, ~, inls] = compute_residuals_l1(A, b, c, d, theta, th);
        disp(['Alpha = ' num2str(alpha) ' P=' num2str(P) ' f=' num2str(F) ' Q=' num2str(Q) ' Inls = ' num2str(numel(inls)) ' maxnorm = ' num2str(max(res1(inls)))]);
        
        disp('---------------------------------------------------------------------');
        if ( Q <= QThresh || alpha > alphaMax)
            break;
        end
        alpha = kappa*alpha;
        
    end
    
    runtime = toc;
    sampleSet = inpSampleSet;
elseif(strcmp(method, 'SS'))   
    [theta, runtime] = SS_Quasi(A, b, c, d, theta0, th, config);
    sampleSet = inpSampleSet;
elseif(strcmp(method, 'IBCO'))   
    [theta, runtime] = IBCO_v2(A, b, c, d, theta0, th, config);
    sampleSet = inpSampleSet;
end
[~, ~, inls]=compute_residuals_l2(A, b, c, d, theta, th,1);

end