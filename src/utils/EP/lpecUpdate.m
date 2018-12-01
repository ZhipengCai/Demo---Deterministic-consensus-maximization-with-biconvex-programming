% Perform Linear Fitting using Linear Equilibrium constraints
%function [theta, inliers, s ,r, iter, currP, f, Q] = lpecUpdate(X, Y, th, theta0, alpha, debug)
function [xn, yn, theta, currP, f, Q, iter, runTime] = lpecUpdate(A, b, x0, y0, alpha, debug)
    %disp([' ----------------Alpha  ' num2str(alpha) ' ------------------ ']);
    if (nargin<4); theta0 = randn(d,1); end;
    if (nargin<5); alpha = 1; end;
    if (nargin<6); debug = 0; end;
    %options=optimoptions('linprog','Algorithm','dual-simplex','Display', 'off');
    % Converting to parameters for MaxFS; |Ax - b|<= th <=>-th<=Ax-b<=th
    %[A, b] = genLinearMaxFSMatrix(X, Y, th);        
    %inls0 = findInliers(theta0);       
    %disp(['Iter=0, Inls  = ' num2str(numel(inls0))]);   
    n = size(A,1);      % Now, the number of constraints is doubled.
    d = size(A,2);           
    maxIter = 1000;        
    % Build the matrices
    %M12 = eye(n); 
    %M21 = -eye(n);
    N11 = [-A A*ones(d,1)];
    %N22 = [];
    %q1 = b;
    %q2 = ones(n,1);
        
    ri = x0(1:n);
    si = x0(n+1:end);
    
    %wi = theta0;
    yi = y0;
    [currP, f, Q] = evalP(ri, si, yi, alpha);
    disp([ 'Init P = ' num2str(currP) ' Init F = ' num2str(f) ' Init Q = ' num2str(Q) ]);
    runTime = tic;
    for iter = 1:maxIter                
        Ays = [-N11 -1*eye(n)];
        Ays = [Ays; -eye(n+d+1)];
        bys = [b; zeros(n+d+1,1)];
        fys = [N11'*ri; ones(n,1)];
        % Stack [w;s] to LP
        %[ys,~,exitFlag]  = linprog(fys, Ays, bys,[],[],[],[],[], options);
        %[ys, exitFlag] = cvxLinProg(fys, Ays, bys);
        %[ys, exitFlag] = lpsolve(fys, Ays, bys);        
        %[ys, exitFlag] = sedumiLinProg(fys, Ays, bys);
        [ys, exitFlag] = gurobiLinProg(fys, Ays, bys);
        %[ys, exitFlag] = lpsLinProg(fys, Ays, bys);
        
        %syy = [ys1 ys]
        if (exitFlag<1)
            disp(exitFlag);
            disp('LP1 failed');          
            break;
        end
        
        yip = ys(1:d+1);
        sip = ys(d+2:end);


        %[P1, F1, Q1] = evalP(ri, sip, yip, alpha);
        %disp (['-----LP 1: P = ' num2str(P1) ' F = ' num2str(F1) ' Q = ' num2str(Q1)]);
        %thetai = computeTheta(yip);
         
         %Test 
         %ss = sip + ri.*(N11*yip);
         %Ab = A*thetai + b;
%         temp = sip - A*thetai + b; temp = temp(1:n);
%         temp2 = sip + A*thetai + b; temp2 = temp2(1:n);
%         temp3 = -sip + A*thetai + b;temp3 = temp3(1:n);
%         temp4 = -sip - A*thetai + b; temp4 = temp4(1:n);
%         At = [sip A*thetai b Ab ];
        
        fr = ones(n,1)+ alpha*(N11*yip+b);
        Ar = [eye(n); -eye(n)];
        br = [ones(n,1); zeros(n,1)];

        %[rip,~,exitFlag] = linprog(fr, Ar, br,[],[],[],[],[], options);
        %[rip, exitFlag] = cvxLinProg(fr, Ar, br);
        %[rip, exitFlag] = lpsolve(fr, Ar, br);        
        %[rip, exitFlag] = sedumiLinProg(fr, Ar, br);
        [rip, exitFlag] = solve_u_lp(fr);       
        %[rip, exitFlag] = lpsLinProg(fr, Ar, br);
        %[rip, exitFlag] = gurobiLinProg(fr, Ar, br);                
        
        [Pnew, Fnew, Qnew] = evalP(rip, sip, yip, alpha);
        
        
          
        if (debug)
            %tti = computeTheta(yi);
            %inlsi = findInliers(tti);
            disp(['Iter = ' num2str(iter) ' P = ' num2str(Pnew) ' F = ' num2str(Fnew) ' Q = ' num2str(Qnew) ' Pdiff =' num2str(Pnew-currP)....
                ]);
        end

        %        ri = rip;
        %if (Pnew>=currP); break; end;      
        if (abs(Pnew-currP)<=1e-9); break; end;

        % Update new results
        ri = rip;
        si = sip;
        yi = yip;


        
        currP = Pnew;
        f = Fnew;
        Q = Qnew;
        
        %if (abs(Q)<=1e-15); break; end;
    end
    
    sc = find(si==0);    
    %theta = yi(1:d) - repmat(yi(d+1),d,1);
    theta = computeTheta(yi);
%     s = si;
%     r = ri;
    xn = [ri; si];
    yn = yi;
    %inliers = findInliers(theta);    
    %disp(theta');
    runTime = toc(runTime);

function [P, F, Q] = evalP(r, s, y, alpha)
    P = sum(abs(r)) + alpha*r'*(s+N11*y+b) + alpha*s'*(-r+ones(n,1));
    F = sum(r);
    %aQ = P - F;
    %Q = aQ/alpha;
    Q = r'*(s+N11*y+b) + s'*(-r+ones(n,1));
end

% function inliers = findInliers(theta)
%     sc1 = consistentSet(theta);
%     n1 = size(X,1);
%     inliers=[];
%     for i=1:length(sc1)
%         if sc1(i)<=n1
%             if (~isempty(find(sc1==sc1(i)+n1)))
%                 inliers = [inliers; sc1(i)];
%             end
%         end
%     end
% end

function cs = consistentSet(theta)
    cs = find(-A*theta+b>=0);
end

function tt = computeTheta(yi)
    tt = yi(1:d) - repmat(yi(d+1),d,1);
end

end


