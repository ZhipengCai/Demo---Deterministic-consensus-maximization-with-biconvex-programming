function [x, flag] = lpsolve(f,A,b)
    
    e = -1*ones(size(b,1),1);
    %lb = -5000*ones(size(b,1),1);
    %up = 5000*ones(size(b,1),1);
    [obj, x] = lp_solve(f,A,b,e);    %1 2 3 4 7 
    
    flag=1;
    if (numel(x)==0); flag = -1; end;
       
    
    %disp(cvx_status);
end


