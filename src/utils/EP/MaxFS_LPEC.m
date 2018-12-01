% Perform Linear Fitting using Linear Equilibrium constraints
function [xn, yn, theta,  P, F, Q, iter] = MaxFS_LPEC(A, b, cd,dd, th, x0, y0, alpha, denLim, lpsolver) %Note: cd, dd are to enforce cheirality contraints

if (nargin < 9)  denLim = 0; end;
if (nargin<8); alpha = 1; end;

n = size(A,1);
d = size(A,2);

% Build the matrices
M12 = eye(n);
M21 = -eye(n);
N11 = [-A A*ones(d,1)];

% generate cN matrix similar to N11 for A to enforce c_i^T\theta + d_i > 0
cdp = cd';
cN11 = [-cdp cdp*ones(d,1)];
cN11 = cN11;
ddp = dd' - denLim*ones(numel(dd),1);

% Return to normal :)
N22 = [];
q1 = b;
q2 = ones(n,1);

options=optimoptions('linprog','Algorithm','dual-simplex','Display', 'off'); %LP Options (only used for MATLAB's linprog)

ri = x0(1:n);
si = x0(n+1:end);
yi = y0;

iter = 0;
[P, F, Q] = evalP(ri, si, yi, alpha);

while (true)
    iter = iter + 1;
    Ays = [-N11 -1*eye(n)];
    Ays = [Ays; -eye(n+d+1)];
    bys = [b; zeros(d+1+n,1)];
    
    %         %Enforcing H33 = 1
    %         Ays = [Ays; [zeros(1,d-1) -1 1 zeros(1,n)]];
    %         bys = [bys; -1];
    %         Ays = [Ays; [zeros(1,d-1) 1 -1 zeros(1,n)]];
    %         bys = [bys; 1];
    
    % Now add the cheirality constraints, only add if user wants denLim > 0
    if (denLim > 0)
        Ays = [Ays; [cN11 zeros(size(cN11,1), n)]];
        bys = [bys; ddp];
        disp(['Adding Constraints CiTx+d > ' num2str(denLim)] );
    end
    
    fys = [N11'*ri; ones(n,1)];
    
    
    % Stack [w;s] to LP
    %[ys,~,exitFlag]  = linprog(fys, Ays, bys,[],[],[],[],[], options);
    %[ys, exitFlag] = cvxLinProg(fys, Ays, bys);
    %[ys, exitFlag] = lpsolve(fys, Ays, bys);
    %[ys, exitFlag] = sedumiLinProg(fys, Ays, bys);
    [ys, exitFlag] = feval(lpsolver, fys, Ays, bys);
    
    
    
    %syy = [ys1 ys]
    if (exitFlag<1)
        disp(exitFlag);
        disp('LP1 failed');
        break;
    end
    
    yip = ys(1:d+1);
    sip = ys(d+2:end);
    
    
    fr = ones(n,1)+ alpha*(N11*yip+b);
    Ar = [eye(n); -eye(n)];
    br = [ones(n,1); zeros(n,1)];
    
    %[rip,~,exitFlag] = linprog(fr, Ar, br,[],[],[],[],[], options);
    %[rip, exitFlag] = cvxLinProg(fr, Ar, br);
    %[rip, exitFlag] = lpsolve(fr, Ar, br);
    %[rip, exitFlag] = sedumiLinProg(fr, Ar, br);
    [rip, exitFlag] = solve_u_lp(fr);
    
    
    if (exitFlag<1)
        disp(exitFlag);
        disp('LP2 failed');
        break;
    end
    
    [Pnew, Fnew, Qnew] = evalP(rip, sip, yip, alpha);
    
    %if (abs(Pnew-P)<=PThresh); break; end;
    if ((Pnew>=P)); break; end;
    % Update new results
    ri = rip;
    si = sip;
    yi = yip;
    
    tti = computeTheta(yi);
    inlsi = consistentSet(tti);
 
    
    P = Pnew;
    F = Fnew;
    Q = Qnew;
    
end

xn = [ri; si];
yn = yi;
theta = computeTheta(yi);


    function [P, F, Q] = evalP(r, s, y, alpha)
        F = sum(r);
        Q = r'*(s+N11*y+b) + s'*(-r+ones(n,1));
        P = F + alpha*Q;
    end


    function cs = consistentSet(theta)
        cs = find(-A*theta+b>=0);
    end

    function tt = computeTheta(yi)
        tt = yi(1:d) - repmat(yi(d+1),d,1);
    end

end


