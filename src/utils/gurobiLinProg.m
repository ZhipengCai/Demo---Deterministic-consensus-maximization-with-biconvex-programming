function [x, flag] = gurobiLinProg(f,A,b)

    n = size(b,1);
    d = size(f,1);
    flag = 1;
    %x = rand(d,1);
    %Construct Gurobi Model
    model.obj = f;
    model.A = sparse(A);
    model.rhs = b;    
    model.sense = '<';
    model.modelsense = 'min';
    %params.outflag = 0;
    params.FeasibilityTol = 1e-9;
    params.IntFeasTol = 1e-9;
    params.OptimalityTol = 1e-9;
    params.MarkowitzTol=1e-4;
    %params.NormAdjust = 0;
    %params.ScaleFlag = 0;
    %params.PerturbValue = 0;
%     params.MIPGapAbs = 0;
%     params.Sifting = 2;
%     params.SiftMethod = 1;
%     params.MIPGap = 0;
%     %params.Quad = 1;
    
    params.OutputFlag = 0;
    params.method  = -1;
    %params.BarIterLimit = 5000;
    results = gurobi(model, params);
    %disp(results);
    x = [];
    if (~(strcmp(results.status,'OPTIMAL') || strcmp(results.status,'SUBOPTIMAL')))
        flag = -1;
    else 
        x = results.x;
    end
    

    
%     c = [f; -f; zeros(n,1)];
%     Ap = [A -A eye(n)];
%     
%     K.l = max(size(Ap));
%     pars.fid = 0;
%     
%     [xxp,~,info] = sedumi(Ap, b, c, K, pars);
%     
%     x = xxp(1:d);
%     xp = xxp(d+1:2*d);
%     x = x-xp;
%     
%     flag = 1;
%     if (info.numerr==2); flag = -1; end;
    
end


