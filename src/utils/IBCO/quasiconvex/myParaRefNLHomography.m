function [theta, runtime] = myParaRefNLHomography(lA, lb, lAA, lbb, c, d, theta0, config)
% Normalize matrix if a 3x3 matrix was supplied
if (length(theta0)==9)
    theta0 = theta0./theta0(end);
    theta0 = theta0(1:end-1);
end

[x0, y0] = genStartingPoint(lA, lb, theta0);
% From theta0, compute u0,s0, v0
%finds the number of inliers
[nInls]=inlierCountQuasiPRNL(lAA, lbb, theta0);
theta = theta0;
NoConstraints = size(lA,1)/4;
NoPoints = numel(d);

tic

%step 0:
NoConstraints
nInls
deltaMax = NoConstraints-nInls;
deltaMin = 0;

[x1,y1, theta0, fDeltaMin] = fwPRNLQuasiconvex(lA, lb, c, d, x0, y0, 0, config);
%Find Inliers and exit if problem solved for 0 outlier situation
if fDeltaMin <= config.QThresh
    theta = theta0;
    runtime = toc;
    return;
end

[nInls]=inlierCountQuasiPRNL(lAA, lbb, theta0);
if deltaMax>(NoConstraints-nInls)
    disp(['noInliers: ' num2str(nInls)]);
    deltaMax = NoConstraints-nInls;
    theta = theta0;
    x0=x1;
    y0=y1;
end

delta = floor(2/3*deltaMax);
while deltaMax>(deltaMin+1)
    %step 1:
    %execute active-set algorithm
    [x1, y1, theta0, fDelta] = fwPRNLQuasiconvex(lA, lb, c, d, x0, y0, delta, config);
    %step 2:
    [nInls]=inlierCountQuasiPRNL(lAA, lbb, theta0);
    deltaMax
    
    if deltaMax>=(NoConstraints-nInls)
        disp(['noInliers: ' num2str(nInls)]);
        deltaMax = NoConstraints-nInls;
        theta = theta0;
        x0=x1;
        y0=y1;
    end
    %step 3:
    if fDelta <= config.QThresh
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
