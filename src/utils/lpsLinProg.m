function [x, flag] = lpsLinProg(f,A,b)

    addpath('../solvers/LPS/');
    n = size(b,1);
    d = size(f,1);
    
    c = [f; -f]';
    Ap = [A -A];
    
%     K.l = max(size(Ap));
%     pars.fid = 0;
    
    [~,xxp,sts, warn] = LP_solve(Ap, b, c, 'min');
    
    x = xxp(1:d);
    xp = xxp(d+1:2*d);
    x = x-xp;
    
    flag = 1;
    if (info.numerr==2); flag = -1; end;
    
end


