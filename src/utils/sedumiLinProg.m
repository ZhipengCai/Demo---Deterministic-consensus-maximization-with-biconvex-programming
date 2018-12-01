function [x, flag] = sedumiLinProg(f,A,b)

    n = size(b,1);
    d = size(f,1);
    
    c = [f; -f; zeros(n,1)];
    Ap = [A -A eye(n)];
    
    K.l = max(size(Ap));
    
    %Setting parameters
    pars.fid = 0;
    pars.maxiter = 500;
    pars.eps = 0;
    
    
    [xxp,~,info] = sedumi(Ap, b, c, K, pars);
    
    x = xxp(1:d);
    xp = xxp(d+1:2*d);
    x = x-xp;
    
    flag = 1;
    if (info.numerr==2); flag = -1; end;
    
end


