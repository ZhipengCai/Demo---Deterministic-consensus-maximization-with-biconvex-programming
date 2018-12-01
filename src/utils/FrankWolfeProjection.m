% Find a projection of a point y to the polyhedron A*x+b>0 using FW
% min ||x-y||^2
% s.t. A*x + b >=0 

function theta = FrankWolfeProjection(y, A, b, x0, maxIter)    
    grad = @(xk) 2*(xk - y);   % compute gradient
    f = @(xk) norm(xk-y)^2;    % compute f(X)        
    x = x0;        
    for iter=0:maxIter            
        grd = grad(x);
        s = gurobiLinProg(grd,-A, b);                       
        if (length(s) < length(x0))            
            theta = x;                    
            return;
        end
         % Perform line serach to find gamma:
        ax = x-y; bx = s - x;
        as = norm(bx)^2; bs = 2*ax'*bx;
        gamma = -bs/(2*as);     
        if (gamma<0 || gamma>1); gamma = 2./(iter+2); end
        prevx = x;              
        dualGap = (x-s)'*grd;
        x = (1-gamma).*x + gamma.*s;
        prevf = f(prevx); fx = f(x);                
        disp([' FW Quadprog iter = ' num2str(iter) ' prev f = ' num2str(prevf) ' f = ' num2str(fx) ' Duality gap = ' num2str(dualGap) ]);                
        if (dualGap < 0.05)
            theta = x; 
            break;            
        end                
    end

end