%frank-wolfe algorithm for parametric formulation
% A, b: input data
% x0, y0: starting point
% lpsolver: either 'sedumi' or 'gurobi'

function [xn, yn, theta, currF] = fwPR(A, b, x0, y0, delta, lpsolver)

    debug  = 1;        % Print debug information during run
    n = size(A,1);     
    d = size(A,2);           
    maxIter = 1e9;     % Safeguard to prevent infinite loop. 
    
    N = [-A A*ones(d,1)];        
    ri = x0(1:n);       % r corresponds to u in the paper
    si = x0(n+1:end);   % s corresponds to s in the paper     
    yi = y0;            % y corresponds to v in the paper
    
    
    disp(['+++Running PR with delta = ' num2str(delta) '......']);
    for iter = 1:maxIter
        % prepare parameteres of LP1
        Ays = [-N -1*eye(n)];
        Ays = [Ays; -eye(n+d+1)];
        bys = [b; zeros(n+d+1,1)];
        fys = [N'*ri; ones(n,1)];        
        
        % LP1: Fix u, solve s and v
        [ys, exitFlag] = feval(lpsolver, fys, Ays, bys);               
        
        if (exitFlag<1)  
            % If for some reason, the solver can't solve the LP completely, 
            % need to check the solver or the data. This rarely happens for normal run
            disp(exitFlag);
            disp('LP1 failed. Please check solver configurations or data');          
            break;
        end
        % Extract s and v from LP1 
        yip = ys(1:d+1);
        sip = ys(d+2:end);
        
        % LP2: Fix s and v, solve u
        fr = (N*yip+b);        
        [rip, ~] = solve_u_lpPR(fr,delta);          % u can be solved in close form        
        [Fnew] = evalF(rip, sip, yip);               
        
        % Update new results
        ri = rip;
        si = sip;
        yi = yip;               
        F = Fnew;
        
        if (debug)
            disp(['Iter = ' num2str(iter) ' F = ' num2str(F) ]);
        end
        
        if (iter>1&&(abs(Fnew-currF)<=1e-9 || Fnew >= currF)); break; end;  % If we reached the stationary point, then break;
        currF = Fnew;
    end
    
    theta = computeTheta(yi);
    xn = [ri; si];
    yn = yi;
    
    function [F] = evalF(r, s, y)
    F = abs(r'*(N*y+b) + s'*ones(n,1));
    end
    
    function tt = computeTheta(yi)
    % From v, compute theta  by substracting the (d+1)-th element from  
    % the first d elements of v
        tt = yi(1:d) - repmat(yi(d+1),d,1);
    end
end