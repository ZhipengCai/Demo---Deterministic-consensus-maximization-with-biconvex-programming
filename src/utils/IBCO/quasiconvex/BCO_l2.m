%biconvex optimization for l2 norm refinement problem
function [sx1,u1, theta, currF,FProj] = BCO_l2(A, b, c, d, modelsx, sx0, u0, delta, config, epsilon, epsilon_ac,n,d1,d2,solver,problem)
if delta == 0 %solve socp when delta = 0
    u1 = zeros(n,1);
    if(strcmp(solver,'gurobi'))
        modelsx.obj = [zeros(1,d2), ones(1,n), zeros(1,d1*n+n)];
        params.outputflag = 0;
        result = gurobi(modelsx, params);
        sx1 = result.x(1:d2+n)
        currF = result.objval
    else
        pars.fid=0;
        bt = [zeros(d2,1); -ones(n,1)];
        currF = -bt'*sx0;
        [~,sx1,info] = sedumi(modelsx.At,bt,modelsx.ct,modelsx.K,pars);
        theta = sx1(1:d2);
        currF = -bt'*sx1;
        %enforce the constraint on theta for fundamental matrix estimation
        if(strcmp(problem,'fun'))
            F = [theta(1:3)';theta(4:6)';[theta(7),1,theta(8)]];
            FPrime = enforceFundamentalConstraint(F);
            theta = [FPrime(1,:) FPrime(2,:) FPrime(3,1) FPrime(3,3)]';
            [~,fr,~] = compute_residuals_l2(A,b,c,d,theta,epsilon);
            k = (d2+1):(d2+n);
            sx1(k) = fr;
            sx1(k(fr<=1e-9)) = 0;
            FProj = -bt'*sx1;
        else
            FProj = currF;
        end
    end
else %use biconvex optimization technique
    maxiter = 1e6;%safe guard to prevent infinite iterations
    currF = 1e9;
    sx1 = sx0;
    k = (d2+1):(d2+n);
    for i = 1:maxiter
        preF = currF;
        %1. fix sx, solve u
        [~,fr,~] = compute_residuals_l2(A,b,c,d,sx1(1:d2),epsilon);
        [u1, ~] = solve_u_lpPR(-fr,delta);          % u can be solved in close form
        %2. fix u, solve sx
        if(strcmp(solver,'gurobi'))
            modelsx.obj = [zeros(1,d2), 1-u1, zeros(1,d1*n+n)];
            params.outputflag = 0;
            result = gurobi(modelsx, params);
            sx1 = result.x(1:d2+n)
            currF = result.objval
        else
            pars.fid=0;
            bt = [zeros(d2,1); u1-1];
            [~,sx1,info] = sedumi(modelsx.At,bt,modelsx.ct,modelsx.K,pars);
            [~,fr,~] = compute_residuals_l2(A,b,c,d,sx1(1:d2),epsilon);
            sx1(k) = fr;
            sx1(k(fr<=1e-9)) = 0;
            currF = sum(sx1(d2+1:d2+n)-u1.*fr);
            theta = sx1(1:d2);
        end
        if(i>1 && (abs(preF-currF)<1e-9 || preF<currF))
            currF = sum(sx1(d2+1:d2+n)-u1.*fr);
            if(strcmp(problem,'fun'))
                F = [theta(1:3)';theta(4:6)';[theta(7),1,theta(8)]];
                FPrime = enforceFundamentalConstraint(F);
                theta = [FPrime(1,:) FPrime(2,:) FPrime(3,1) FPrime(3,3)]';
                [~,fr,~] = compute_residuals_l2(A,b,c,d,theta,epsilon);
                sx1(k) = fr;
                sx1(k(fr<=1e-9)) = 0;
                FProj = sum(sx1(d2+1:d2+n)-u1.*fr);
            else
                FProj = currF;
            end
            break;
        end
    end
end
end