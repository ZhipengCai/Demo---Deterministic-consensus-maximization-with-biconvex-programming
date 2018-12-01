% solve quadratic problem with one constraint using bisection
% min ||x - z||_2^2
% s.t. x^T*Aq*x + bq^T * x + r <= 0
% Refer to Appendix B (http://stanford.edu/~boyd/papers/pdf/qcqp.pdf)
% This code was converted from python from
% (https://github.com/cvxgrp/qcqp/blob/master/qcqp/utilities.py)

function [theta, cost] = qcqp_oneconstr(z, Aq, bq, r, Q, D, tol)
   
    if (nargin<7); tol = 1e-3; end
           
    % If z satisfies the constraint, it is the solution
    if z'*Aq*z + bq'*z + r <= 0
        theta = z; cost = 0;
        return
    end
        
     % Eigen Decomposition for Aq:
%     P = Aq;   %q = bq; % Notation used in paper
%     Psymm = 0.5*(P + P');
%     [Q,D] = eig(Psymm);   
    lambda = D*ones(size(D,1),1);
     %I = eye(size(D,1));
    
    q = bq;
    zhat = Q'*z;
    qhat = Q'*q;
    
    % Solve the problem with the transformed variables 
    % min_xhat ||xhat - zhat||_2^2
    % s.t. sum(lamdba_i xhat_i^2) + qhat^T xhat + r = 0 
    
    %xhat = @(nu) -0.5.*inv((I + nu.*D))*(nu*qhat - 2*zhat);
    %xhat2 = @(nu) -0.5.*(nu*qhat - 2*zhat)./((I + nu.*D)*ones(size(D,1),1));
    xhat = @(nu) -0.5.*(nu*qhat - 2*zhat)./(nu.*lambda+ones(length(lambda),1));
    phi = @(xh) xh'*D*xh + qhat'*xh + r;
    
    s = -inf;
    e = inf;
%     
     for i=1:length(lambda)
         if (lambda(i)> 0); s = max(s, -1./lambda(i)); end
         if (lambda(i)< 0); e = min(e, -1./lambda(i)); end
     end
%     
      
     if (s == -inf)
         s = -1.0;
         while phi(xhat(s))<=0; s = s*2.0; end        
     end
%     
     if (e == inf)
         e = 1.0;
         while phi(xhat(e))>=0; e = e*2.0; end
     end
    %plotPhi();
    while (e-s > tol)
        m = (s+e)/2.0;
        p = phi(xhat(m));        
        if abs(p)<tol % m is the solution            
             s = m;
             e = m;          
             break;        
        elseif p > 0;  s = m;
        elseif p < 0; e = m;
        else 
            s = m;
            e = m;          
            break;        
        end                                     
    end
    nu = 0.5*(s+e);
    xh = xhat(nu);
    theta = Q*xh;      
    cost = norm(theta-z,2)^2;
%     tt = theta;
%     t1 = tt'*Aq*tt + bq'*tt + r
%     disp('---');
    
%     function [Q, D, lambda, s, e] = bisection_param(Aq)
%     [Q,D] = eig(Aq);
%     lambda = D*ones(size(D,1),1);
%     s = -inf;
%     e = inf;    
%     for i=1:length(lambda)
%         if (lambda(i)> 0); s = max(s, -1./lambda(i)); end
%         if (lambda(i)< 0); e = min(e, -1./lambda(i)); end
%     end
%     if (s == -inf)
%         s = -1.0;
%         while phi(xhat(s))<=0; s = s*2.0; end        
%     end
%     
%     if (e == inf)
%         e = 1.0;
%         while phi(xhat(e))>=0; e = e*2.0; end
%     end
%     end

    % Test    
%     function plotPhi()
%         close all;
%         nn = s+0.001:0.00001:e-0.001; pp = zeros(size(nn));
%         for i=1:length(nn)
%             pp(i) = phi(xhat(nn(i)));
%         end
%         zeroline = 0.*nn;
%         plot(nn, pp, '.-');
%         hold on;
%         plot(nn, zeroline, '--');
%         
%         
%     end
    
end