function [resn, s, inliers] = compute_residuals_l1(A, b, c, d, x, th,  normalizeH)

    if (nargin < 7); normalizeH = 1; end;        
    if (nargin<6);   th = 1;   end;    
    
    if (normalizeH == 1) && (length(x)==9); x = x./x(end); end; % For homography: If normalize, normalize homography matrix to set H33 = 1;
    
    if (size(A,2)<length(x))
        x=x(1:size(A,2));
    end
    
    nbrimages = numel(d); 
    p = size(A, 2); 
    AxPb = reshape(A*x,2,nbrimages) + b; % Ax + b
    Sqrt_AxPb = (sum(abs(AxPb)));                     % sqrt(Ax + b)
    CxPd = (x'*c + d);
    
    resn = zeros(nbrimages, 1);
    id = abs(CxPd)>0.01;    
    resn(id) = Sqrt_AxPb(id)./abs(CxPd(id));
    resn(~id) = 100000000*max(resn);
    inliers = find(abs(resn) <= th); 
    s = (Sqrt_AxPb./th - CxPd)';
end
