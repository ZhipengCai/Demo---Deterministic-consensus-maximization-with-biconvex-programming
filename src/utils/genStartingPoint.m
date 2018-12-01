function [x0, y0] = genStartingPoint(A, b, theta0)

    n = size(A,1);
    d = size(A,2);
    % Get the initial variables (x0, y0)
    r = zeros(n,1);
    s = zeros(n,1);
    res = -A*theta0 + b;   % Compute residual. If all are consistent, all must be non-negative.
    
    failedIdx = find(res<0);
    r(failedIdx)=1.0;           
    
    s(failedIdx) = -1.0*res(failedIdx);
    
    psi = min(theta0);
    if (psi>=0); psi = 0; else psi=-psi; end;

    % Starting point:
    x0 = [r; s];
    y0 = [theta0+repmat(psi, d, 1);psi];

end