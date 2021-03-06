%reweighted l1 methods (IRLP)
function [theta, inls, nInls, runtime] = SS(x, y, th, theta0, config)
tic
%initialization
theta = theta0;
s = max(0,abs(y-x*theta0)-th);

gamma = config.gammaSS;


w = 1./(gamma+s);
while true
    %solve one LP
    fs = [zeros(length(theta),1);w]; %obj
    As = [x; -x; zeros(size(x))];
    Is = [-eye(numel(y)); -eye(numel(y)); -eye(numel(y))];
    As = [As, Is];
    bs = [th+y; th-y; zeros(numel(y),1)];
    s0=s;
    if isequal(config.solver.LP,@linprog)
        options = optimset('linprog');
        options.Display = 'off';
        [ys, ~, exitFlag] = feval(config.solver.LP, fs,As, bs,[],[],[],[],options);
    else
        [ys, exitflag] = feval(config.solver.LP,fs,As, bs);
    end
    theta = ys(1:length(theta));
    
    s = ys(length(theta)+1:end);
    dif = norm(s-s0);
    w = 1./(gamma+s);
    if(dif<config.QThresh)
        break;
    end
end
inls = [];
nInls = [];
runtime = toc;
end
