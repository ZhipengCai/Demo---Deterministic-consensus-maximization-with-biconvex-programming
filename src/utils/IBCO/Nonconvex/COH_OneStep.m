%COH step for a fixed delta
function [sx1, u1, theta, currF, FProj] = COH_OneStep(A,b,c,d,modelYamip,sx0,u0,delta,config,th,th_accu,n,d1,d2,constraintType)
disp(['COH with delta = ', num2str(delta)]);
if delta == 0 %solve sdp when delta = 0
    %give every point an inlier flag
    u1 = zeros(n,1);
    %setting objective
    Objective = sum(modelYamip.s);
    %yamip solver options
    options = sdpsettings('solver','sedumi','verbose',0, 'allownonconvex', 0);
    sol = optimize(modelYamip.Constraints,Objective,options);
    %extract solution
    if sol.problem == 0
        theta = value(modelYamip.x);
        s = value(modelYamip.s);
        currF = sum(s);
        sx1 = [theta;s];
    else
        display('fail to solve sdp');
        sol.info
        yalmiperror(sol.problem)
    end
    %enforce the constraint
    if strcmp(constraintType,'so(3)')
        %enforce the rotation constraint
        theta = enforceSO3Constraint(theta);
        [resn,s,inls] = compute_residuals_l2(A',b,c,d,theta,th);
        s(s<=config.QThresh) = 0;
        sx1 = [theta;s];
        FProj = sum(s);
    else
        FProj = currF;
    end
    disp(['currF = ', num2str(currF), '; FProj = ',num2str(FProj)]);
else
    maxIter = 1e3;%safe guard to prevent infinite iterations
    sx1 = sx0;
    
    options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0);
    
    %1. fix sx, solve u
    [~,fr,~] = compute_residuals_l2(A',b,c,d,sx1(1:d1),th);
    [u1, ~] = solve_u_lpPR(-fr,delta);
    
    currF = (1-u1)'*sx1(d1+1:d1+n);

    for i = 1:maxIter
        preF = currF;
        
        %2. fix u, solve sx
        %compute current inliers
        Objective = (1-u1)'*modelYamip.s;
        sol = optimize(modelYamip.Constraints,Objective,options);
        %extract solution
        if sol.problem == 0
            theta = value(modelYamip.x);
            s = value(modelYamip.s);
            
            %1. fix sx, solve u
            [~,fr,~] = compute_residuals_l2(A',b,c,d,theta,th);
            [u1, ~] = solve_u_lpPR(-fr,delta);
            
            currF = (1-u1)'*s;
            sx1 = [theta;s];
        else
            display('fail to solve sdp');
            sol.info
            yalmiperror(sol.problem)
        end
        %converge criteria
        if(i>1 && (abs(preF-currF)<1e-9 || preF<currF))
            %enforce the constraint
            if strcmp(constraintType,'so(3)')
                %enforce the rotation constraint
                theta = enforceSO3Constraint(theta);
                %fix sx and solve u at last if enforcing constraint
                [resn,s,inls] = compute_residuals_l2(A',b,c,d,theta,th);
                [u1, ~] = solve_u_lpPR(-s,delta);
                s(s<=config.QThresh) = 0;
                sx1 = [theta;s];
                FProj = (1-u1)'*s;
            else
                FProj = currF;
            end
            disp(['converge at iter: ', num2str(i),'; currF = ', num2str(currF), '; FProj = ',num2str(FProj)]);
            break;
        end
    end
end
end
