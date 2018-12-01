function [fDeltaVec, runtimeVec,  fDeltaVec3, runtimeVec3] = myParaRefTest(x, y, th, theta0, config)
     	[A, b]= genLinearMaxFSMatrix(x, y, th-config.QThresh);         % Collect data into matrix A and vector b
        % to form set of constrains Ax <= b
        [x0, y0] = genStartingPoint(A, b, theta0); % From theta0, compute u0, s0, v0
        
        %finds the number of inliers
        NoPoints = length(y);
        inls = find(abs(y-x*theta0)<=th);
        nInls = numel(inls);
        theta = theta0;

        %step 0:
        deltaMax = NoPoints-nInls;
        deltaMin = 0;

        [x1, y1, theta0, fDeltaMin] = BiAS_StBl(A, b, x0, y0, 0, config.lpsolver);
        %Find Inliers and exit if problem solved for 0 outlier situation
        if fDeltaMin <= config.QThresh
            theta = theta0;
            nInls = NoPoints;
            inls = find(abs(y-x*theta0)<=th);
            return;
        end
	%judge whether current solution is better
        inls = find(abs(y-x*theta0)<=th);
        nInls = numel(inls);
        %step 2:
        if deltaMax>(NoPoints-nInls)
            deltaMax = NoPoints-nInls;
            disp(['noInliers ' num2str(nInls)]);
            theta = theta0;
            x0=x1;
            y0=y1;
        end
	%initialize delta
        delta = floor(2/3*deltaMax);
        tic
        runtimeVec = [0];
        [x2, y2, theta2] = BiAS_StBl_M(A, b, x1, y1, delta, config.lpsolver);
        runtimeVec = [runtimeVec toc];
        
        tic
        runtimeVec3 = [0];
        [x3, y3, theta3] = BiAS_fwBlTest(A, b, x1, y1, delta, config.lpsolver);
        runtimeVec3 = [runtimeVec3 toc];
        
        fDeltaVec3 = [];
        fDeltaVec = [];
        %[x1, y1, theta0, fDeltaVec, runtimeVec] = BiAS_StBlTest(A, b, x0, y0, delta, config.lpsolver);
        %[x2, y2, theta2, fDeltaVec2, runtimeVec2] = BiSimplex_StBlTest(A, b, x0, y0, delta, config.lpsolver);
        %[x3, y3, theta3, fDeltaVec3, runtimeVec3] = BiAS_fwBlTest2(A, b, x0, y0, delta, config.lpsolver); 
end
