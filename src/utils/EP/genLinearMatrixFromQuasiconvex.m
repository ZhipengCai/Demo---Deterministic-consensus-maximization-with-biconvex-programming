% This function generate convert a quasiconvex problem (in L1/Linf) into a
% linear fitting problem
function [lA, lb] = genLinearMatrixFromQuasiconvex(A, b, c, d, th)

  denLim = 0.5;
  thetaNormLim = 1.59;
  lA = [];
  lb = [];
  
  for i=1:2:(size(A,1)-1)  
    a1 = A(i,:);
    a2 = A(i+1,:);  
    idx = ceil(i/2);
    b1 = b(1,idx);
    b2 = b(2,idx);
    lA = [lA;  a1+a2 - th*c(:,idx)']; lb = [lb; th*d(idx) - b1 - b2];
    lA = [lA;  a2-a1 - th*c(:,idx)']; lb = [lb; th*d(idx) - b2 + b1]; 
    lA = [lA;  a1-a2 - th*c(:,idx)']; lb = [lb; th*d(idx) - b1 + b2];
    lA = [lA; -a1-a2 - th*c(:,idx)']; lb = [lb; th*d(idx) + b1 + b2];           
  end      
end


