function [theta, inls, nInls, runtime] = myParaRef_ibcoType(x, y, th, theta0, config)
     	[A, b]= genLinearMaxFSMatrix(x, y, th-config.QThresh);         % Collect data into matrix A and vector b
        % to form set of constrains Ax <= b
        [x0, y0] = genStartingPoint(A, b, theta0); % From theta0, compute u0, s0, v0
        
        %finds the number of inliers
        NoPoints = length(y);
        inls = find(abs(y-x*theta0)<=th);
        nInls = numel(inls);
        theta = theta0;
        tic
        %step 0:
        deltaMax = NoPoints-nInls;
        deltaMin = 0;

        [x1, y1, theta0, fDeltaMin] = fwPR_improved(A, b, x0, y0, 0, config.lpsolver);
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

        while deltaMax>(deltaMin+1)
            %step 1:
            [x1, y1, theta0, fDelta] = fwPR_improved(A, b, x0, y0, delta, config.lpsolver);
            inls = find(abs(y-x*theta0)<=th);
            nInls = numel(inls);
            %step 2:
            if deltaMax>=(NoPoints-nInls)
                deltaMax = NoPoints-nInls;
                disp(['noInliers ' num2str(nInls)]);
                theta = theta0;
                x0=x1;
                y0=y1;
            end
            
            if fDelta <= config.QThresh
                delta = floor((deltaMin+deltaMax)/2);
            else
                p = floor(delta - fDelta*(delta-deltaMin)/(fDelta-fDeltaMin));
                deltaMin = delta;
                fDeltaMin = fDelta;
                if p>deltaMin&&p<deltaMax
                    delta = p;
                else
                    %disp('delta exceeds the range');
                    delta = floor((deltaMin+deltaMax)/2);
                end
            end
        end
        %inls = find(abs(y-x*theta)<=th);
        %nInls = numel(inls);
        runtime = toc;
end
