%frank-wolfe algorithm for homography data using parametric reformulation method 
function [xn, yn, theta,  currF] = fwPRIMQuasiconvex(A, b, cd,dd, x0, y0, delta, config, denThresh) %Note: cd, dd are to enforce strictly cheirality contraints

    if (nargin < 9)  denThresh = 0; end;  
    
    n = size(A,1);
    d = size(A,2);
    % Build the matrices
    
    N11 = [-A A*ones(d,1)];
    % generate cN matrix similar to N11 for A to enforce c_i^T\theta + d_i
    % > 0 to prevent possible divide by 0; 
    cdp = cd';    
    cN11 = [-cdp cdp*ones(d,1)];    
    ddp = dd' - denThresh*ones(numel(dd),1);        
    ri = x0(1:n);
    si = x0(n+1:end);
    yi = y0;            
    iter = 0;

    disp(['+++Running PR with delta = ' num2str(delta) '......']);
    while (true)
        iter = iter + 1; 
        Ays = [-N11 -1*eye(n)];
        Ays = [Ays; -eye(n+d+1)];
        bys = [b; zeros(d+1+n,1)];  
        
        % Now add the cheirality constraints, only add if user wants denThresh > 0
        if (denThresh > 0)
            Ays = [Ays; [cN11 zeros(size(cN11,1), n)]];
            bys = [bys; ddp];
            disp(['Adding Constraints CiTx+d > ' num2str(denThresh)] );
        end   
        
        fys = [zeros(d+1,1); 1-ri];   
        %LP1: fix u, solve s and v
        [ys, exitFlag] = feval(config.lpsolver, fys, Ays, bys);      
        
        %syy = [ys1 ys]
        if (exitFlag<1)
            disp(exitFlag);
            disp('LP1 failed');          
            break;
        end
        %extract s and v from LP1
        yip = ys(1:d+1);
        sip = ys(d+2:end);       
        
        %LP2: Fix s and v, solve u
        fr = N11*yip+b;        
        [rip, ~] = solve_u_lpPR_Homo(fr,delta);                                
        [Fnew] = evalF(rip, sip, yip);
        
        % Update new results
        ri = rip;
        sumUi = sum(ri)
        si = sip;
        yi = yip;               
    
        disp(['Iter = ' num2str(iter) ' F = ' num2str(Fnew) ]);               
        if ( iter>1&&((abs(Fnew-currF)<1e-9) || (Fnew>=currF))); break; end;
        currF = Fnew;                       
    end    
        
    xn = [ri; si];
    yn = yi;
    theta = computeTheta(yi);
    

    function [F] = evalF(r, s, y)        
        F = abs(r'*(N11*y+b) + s'*ones(n,1));      
    end

    function tt = computeTheta(yi)
        tt = yi(1:d) - repmat(yi(d+1),d,1);
    end

end


