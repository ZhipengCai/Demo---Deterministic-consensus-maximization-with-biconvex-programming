function [H, inliers, trials, T1, T2, sampleSet] = ransacfithomography(x1, x2, t, normalize, updateMethod, inpSampleSet, normType, loransacInlThreshold, matchingScores)

    if (nargin < 8)  loransacInlThreshold = 0; end
    if (nargin < 9)  matchingScores = []; end
    if (nargin < 4)  normalize = true;   end
    if (nargin < 7)  normType = 'l2'; end
    
    if ~all(size(x1)==size(x2))
        error('Data sets x1 and x2 must have the same dimension');
    end
    
    [rows,npts] = size(x1);
    if rows~=2 & rows~=3
        error('x1 and x2 must have 2 or 3 rows');
    end
    
    if npts < 4
        error('Must have at least 4 points to fit homography');
    end
    
    T1 = [];  T2 = [];
    if rows == 2    % Pad data with homogeneous scale factor of 1
        x1 = [x1; ones(1,npts)];
        x2 = [x2; ones(1,npts)];        
    end
    
    s = 4;  % Minimum No of points needed to fit a homography.
    
    fittingfn = @homography2d;
    if(strcmp(normType,'l1'))
        distfn    = @homogdist2d_l1;
    else
        distfn    = @homogdist2d_l2;    
    end
    degenfn   = @isdegenerate;
    feedback = 0;
    % x1 and x2 are 'stacked' to create a 6xN array for ransac
    [H, inliers, trials, sampleSet] = myHomographyRansac([x1; x2], fittingfn, distfn, degenfn, s, t, updateMethod, inpSampleSet, feedback, 1000, 10000000000, loransacInlThreshold, matchingScores);
    
    
    theta0 = H(:);
    [A, b, c, d] = genMatrixHomography(x1, x2);
    [~, ~, inls0]=compute_residuals_l2(A, b, c, d, theta0, t,1);
    x = [x1;x2];
    H = homography2d(x(:,inls0));
    H = H';
    theta = H(:);
    [~, ~, inliers]=compute_residuals_l2(A, b, c, d, theta, t,1);
    if(numel(inls0)>numel(inliers))
        theta = theta0;
        inliers = inls0;
    end
    
    
%----------------------------------------------------------------------
% Function to evaluate the l2 symmetric transfer error of a homography with
% respect to a set of matched points as needed by RANSAC.

function [inliers, H] = homogdist2d_l2(H, x, t)
    
    x1 = x(1:3,:);   % Extract x1 and x2 from x
    x2 = x(4:6,:);        
    
    [A, b, c, d] = genMatrixHomography(x1, x2);
    
    H = H'; 
    [~, ~, inliers] = compute_residuals_l2(A, b, c, d, H(:), t,1);    
    
    
%----------------------------------------------------------------------
% Function to determine if a set of 4 pairs of matched  points give rise
% to a degeneracy in the calculation of a homography as needed by RANSAC.
% This involves testing whether any 3 of the 4 points in each set is
% colinear. 
     
function r = isdegenerate(x)

    x1 = x(1:3,:);    % Extract x1 and x2 from x
    x2 = x(4:6,:);    
    
    r = ...
    iscolinear(x1(:,1),x1(:,2),x1(:,3)) | ...
    iscolinear(x1(:,1),x1(:,2),x1(:,4)) | ...
    iscolinear(x1(:,1),x1(:,3),x1(:,4)) | ...
    iscolinear(x1(:,2),x1(:,3),x1(:,4)) | ...
    iscolinear(x2(:,1),x2(:,2),x2(:,3)) | ...
    iscolinear(x2(:,1),x2(:,2),x2(:,4)) | ...
    iscolinear(x2(:,1),x2(:,3),x2(:,4)) | ...
    iscolinear(x2(:,2),x2(:,3),x2(:,4));
    









