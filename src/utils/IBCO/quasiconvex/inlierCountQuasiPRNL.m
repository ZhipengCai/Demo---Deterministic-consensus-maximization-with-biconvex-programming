% Compute resiuals and inlier for homography problem
function [inliers] = inlierCountQuasiPRNL(lA, lb, x, normalizeH)

    if (nargin < 4); normalizeH = 1; end;              
    
    if (normalizeH == 1) && (length(x)==9); x = x./x(end); end; % For homography: If normalize, normalize homography matrix to set H33 = 1;
    
    if (size(lA,2)<length(x))
        x=x(1:size(lA,2));
    end
    
    inls = lA*x<=lb;
    sizeI = length(lb)/4;
    i = 1:sizeI;
    inlsI = inls(i*4)+inls(i*4-1)+inls(i*4-2)+inls(i*4-3);
    inliers = sum(inlsI<4);
end
