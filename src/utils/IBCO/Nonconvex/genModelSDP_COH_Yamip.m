function [modelYamip] = genModelSDP_COH_Yamip(A,b,c,d,th,n,d1,d2,constraintType)
modelYamip.x = sdpvar(d1,1);
modelYamip.s = sdpvar(n,1);
%adding the relaxed non-convex constraints
if strcmp(constraintType,'so(3)') %rotation
    L = eye(4)+conv_SO3(modelYamip.x,'col');
    modelYamip.Constraints = [L>=0];
else
    modelYamip.Constraints = [];
end
%s_i>=0
for i = 1:n
    A_i = A(:,d2*(i-1)+1:d2*i);
    b_i = double(b(:,i));
    if(numel(d) == 0)
        modelYamip.Constraints = [modelYamip.Constraints, ...
                                  modelYamip.s(i)>=0,...
                                  modelYamip.s(i)-norm(A_i'*modelYamip.x+b_i)/th+1>=0];
    else
        c_i = c(:,i);
        d_i = d(i);
        modelYamip.Constraints = [modelYamip.Constraints,...
                                  modelYamip.s(i)>=0,...
                                  modelYamip.s(i)-norm(A_i'*modelYamip.x+b_i)/th+(c_i'*modelYamip.x+d_i)>=0];
    end
end
end