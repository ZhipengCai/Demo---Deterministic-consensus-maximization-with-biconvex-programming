% Frank-Wolfe algorithm
% A, b: input data
% x0, y0: starting point
% alpha: pealty parameter
% lpsolver: either 'sedumi' or 'gurobi'

function [xn, yn, theta, currP, F, Q] = fw(A, b, x0, y0, alpha, lpsolver)

debug  = 1;        % Print debug information during run
n = size(A,1);
d = size(A,2);
maxIter = 1e9;     % Safeguard to prevent infinite loop.

N = [-A A*ones(d,1)];
ri = x0(1:n);       % r corresponds to u in the paper
si = x0(n+1:end);   % s corresponds to s in the paper
yi = y0;            % y corresponds to v in the paper
[currP, F, Q] = evalP(ri, si, yi, alpha);    % Compute the curent P, F, Q

disp(['+++Running EP with alpha = ' num2str(alpha) '......']);
for iter = 1:maxIter
    % prepare parameteres of LP1
    Ays = [-N -1*eye(n)];
    Ays = [Ays; -eye(n+d+1)];
    bys = [b; zeros(n+d+1,1)];
    fys = [N'*ri; ones(n,1)];
    
    % LP1: Fix u, solve s and v
    if isequal(lpsolver,@linprog)
        options = optimset('linprog');
        options.Display = 'off';
        [ys, ~,exitFlag] = feval(lpsolver, fys, Ays, bys,[],[],[],[],options);
    else
        [ys, exitFlag] = feval(lpsolver, fys, Ays, bys);
    end
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
    fr = ones(n,1)+ alpha*(N*yip+b);
    [rip, ~] = solve_u_lp(fr);          % u can be solved in close form
    [Pnew, Fnew, Qnew] = evalP(rip, sip, yip, alpha);
    
    % Update new results
    ri = rip;
    si = sip;
    yi = yip;
    F = Fnew;
    Q = Qnew;
    
    if (debug)
        disp(['Iter = ' num2str(iter) ' P = ' num2str(Pnew) ' F = ' num2str(F) ' Q = ' num2str(Q) ]);
    end
    
    if (abs(Pnew-currP)<=1e-5 || Pnew >= currP); break; end;  % If we reached the stationary point, then break;
    currP = Pnew;
end

theta = computeTheta(yi);
xn = [ri; si];
yn = yi;


    function [P, F, Q] = evalP(r, s, y, alpha)
        P = sum(abs(r)) + alpha*r'*(s+N*y+b) + alpha*s'*(-r+ones(n,1));
        F = sum(r);
        Q = abs(r'*(s+N*y+b) + s'*(-r+ones(n,1)));
    end


    function tt = computeTheta(yi)
        % From v, compute theta  by substracting the (d+1)-th element from
        % the first d elements of v
        tt = yi(1:d) - repmat(yi(d+1),d,1);
    end

end


