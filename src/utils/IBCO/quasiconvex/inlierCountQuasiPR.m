% Compute number of inlier (number of satisfied constraints for l1 norm) for homography problem
function [inliers] = inlierCountQuasiPR(lA, lb, x, normalizeH)

    if (nargin < 4); normalizeH = 1; end;              
    
    if (normalizeH == 1) && (length(x)==9); x = x./x(end); end; % For homography: If normalize, normalize homography matrix to set H33 = 1;
    
    if (size(lA,2)<length(x))
        x=x(1:size(lA,2));
    end
    

    inliers = sum(lA*x<=lb);
end
