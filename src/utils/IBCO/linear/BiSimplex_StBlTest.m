%active set algorithm using heuristics on ui and steepest edge pivoting
%rule and Bland's rule for anti-cycling
function [xn, yn, theta, currFVec, runtimeVec] = BiSimplex_StBlTest(A, b, x0, y0, delta, lpsolver)
%if delta = 0 solve the lp produced by setting all ui = 0
if delta == 0
    n = size(A,1);
    d = size(A,2);
    m = d+1;
    
    N = [-A A*ones(d,1)];
    ri = zeros(n,1);       % r corresponds to u in the paper
    si = x0(n+1:end);   % s corresponds to s in the paper
    yi = y0;            % y corresponds to v in the paper
    
    Ays = [-N -1*eye(n)];
    Ays = [Ays; -eye(n+d+1)];
    bys = [b; zeros(n+d+1,1)];
    fys = [zeros(d+1,1); (1-ri)];
    
    % LP1: Fix u, solve s and v
    %[ys, ~, exitFlag] = feval(lpsolver, fys, Ays, bys);
    [ys, exitFlag] = feval(lpsolver, fys, Ays, bys);
    if (exitFlag<1)
        % If for some reason, the solver can't solve the LP completely,
        % need to check the solver or the data. This rarely happens for normal run
        disp(exitFlag);
        disp('LP1 failed. Please check solver configurations or data');
    end
    % Extract s and v from LP1
    yi = ys(1:d+1);
    si = ys(d+2:end);
else
    n = size(A,1);
    d = size(A,2);
    m = d+1;
    maxIter = 1e9;     % Safeguard to prevent infinite loop.
    accuracy = 1e-9;
    
    N = [-A A*ones(d,1)];
    yi = y0;            % y corresponds to v in the paper
    
    disp(['+++Running PR_VE with delta = ' num2str(delta) '......']);
    %give an initial vertex solution using the starting point v0
    %run one iteration of improved frank-wolfe(think about other ways if possible)
    % prepare parameteres of LP1
    % Extract s and v from LP1
    
    % LP2: Fix s and v, solve u
    fr = (N*yi+b);
    [ri, ~] = solve_u_lpPR(fr,delta);          % u can be solved in close form
    
    
    Ays = [-N -1*eye(n)];
    Ays = [Ays; -eye(n+d+1)];
    bys = [b; zeros(n+d+1,1)];
    fys = [zeros(d+1,1); (1-ri)];
    
    % LP1: Fix u, solve s and v
    %[ys, ~,exitFlag] = feval(lpsolver, fys, Ays, bys);
    [ys, exitFlag] = feval(lpsolver, fys, Ays, bys);
    if (exitFlag<1)
        % If for some reason, the solver can't solve the LP completely,
        % need to check the solver or the data. This rarely happens for normal run
        disp(exitFlag);
        disp('LP1 failed. Please check solver configurations or data');
    end
    si = ys(d+2:end);
    yi = ys(1:d+1);
    
    %start the iteration
    currF = evalF(ri,si,yi)
    if(currF<1e-9)
        %compute final result
        theta = computeTheta(yi);
        xn = [ri; si];
        yn = yi;
        return;
    end
    
    %initialize active set and other parameters for set U & (S,V)
    %current constraints
    AX = sparse([eye(n), N; eye(m+n)]');
    bX = sparse([-b; zeros(m+n,1)]);
    
    Xi = [si;yi];
    YiOne = find(ri==1);
    YiZero = find(ri==0);
    
    %current derivative
    derFX = [ones(n,1); N'*ri];
    derFY = b+N*yi;
    
    [valLeaveY, idxLeaveY] = min(derFY(YiZero));
    [valEnterY, idxEnterY] = min(-derFY(YiOne));
    %number of current outliers
    noOut = length(YiOne);
    
    
    if noOut == delta && valLeaveY+valEnterY<valEnterY
        valReduceY = valLeaveY+valEnterY;
        changeModeY = 0;
    elseif noOut == delta
        valReduceY = valEnterY;
        changeModeY = 2;
    elseif valLeaveY<valEnterY
        valReduceY = valLeaveY;
        changeModeY = 1;
    else
        valReduceY = valEnterY;
        changeModeY = 2;
    end
    
    %active set(elements with value 1 represents the active set)
    aIdxX = ((AX'*Xi - bX)<accuracy);
    %index to every active constraint
    aIndexX = find(aIdxX==1);
    %should check for degeneracy: if degeneracy happens, use bland's rule instead
    if length(aIndexX)>m+n
        disp('degeneracy happens for s and v')
        bland = 1;
    else
        bland = 0;
    end
    %matrix for active set
    aAX = AX(:,aIdxX);
    if(bland == 1)
        %find linear independent submatrix
        [aAX,aSubIndexX] = lisubmat(aAX,'cols');
        aIndexX = aIndexX(aSubIndexX);
        aIdxX = zeros(length(aIdxX),1);
        aIdxX(aIndexX) = 1;
    end
    %abX = bX(aIdxX);
    
    %inverse matrix for active set
    aInvAX = inv(aAX);
    %instead of computing inverse, keep LU and solve linear equations for
    %lamda and s
    %     [LX,UX,PX] = lu(aAX);
    
    %flag showing that s,v update has converged
    convergeX = 1;
    %calculate whether there is possible improvements on ui
    if valReduceY<0
        convergeY = 0;
    else
        convergeY = 1;
    end
    %     opts1.LT = true, opts1.TRANSA = false;
    %     opts2.UT = true, opts2.TRANSA = false;
    
    %start testing
    currF = evalF(ri,si,yi)
    tic;
    currFVec = [currF];
    runtimeVec = [toc];
    
    for iter = 1:maxIter
        %compute multiplier
        lambdaX = aInvAX*derFX;
        %         %instead of using inverse matrix, use LU and solve linear equations
        %         lambdaXBar = linsolve(LX,derFX,opts1);
        %         lambdaX = linsolve(UX,lambdaXBar,opts2);
        
        
        %find all the negative lambda
        idxNegLambdaX = find(lambdaX<-accuracy);
        %judge whether the algorithm has converged
        if isempty(idxNegLambdaX)
            convergeX = 1;
            leaveLambdaX = inf;
        else
            %calculate the squared norm of aInvAX for each negative lambda
            if bland == 0
                negAInvAX = aInvAX(idxNegLambdaX,:);
                sqNorm = sum(negAInvAX.*negAInvAX,2);
                [convSpeed, idxMinNegLambdaX] = max((lambdaX(idxNegLambdaX).*lambdaX(idxNegLambdaX))./sqNorm);
                idxMinLambdaX = idxNegLambdaX(idxMinNegLambdaX);
                leaveLambdaX = lambdaX(idxMinLambdaX);
            else
                idxMinNegLambdaX = 1;
                idxMinLambdaX = idxNegLambdaX(idxMinNegLambdaX);
                leaveLambdaX = lambdaX(idxMinLambdaX);
                convSpeed = 0;
            end
        end
        
        
        if(convergeX==1 && convergeY==1)
            disp('converged');
            iter
            break;
        end
        
        convSpeedY = valReduceY*valReduceY;
        if(changeModeY == 0)
            convSpeedY = convSpeedY*0.5;
        end
        
        if convergeX == 0 && (convergeY==1||convSpeed<convSpeedY)
            %disp('entering x');
            if bland == 0
                [aIdxX, aIndexX, Xi, aInvAX, bland] = pivotVE_Deg(AX,bX,Xi,aIdxX,aIndexX,aInvAX,idxMinLambdaX);
            else
                [aIdxX, aIndexX, Xi, aInvAX, convergeX, bland] = pivotVE_Bland(AX,bX,Xi,aIdxX,aIndexX,aInvAX,derFX,idxNegLambdaX);
            end
            si = Xi(1:n);
            yi = Xi((n+1):end);
            derFY = b+N*yi;
            [valLeaveY, idxLeaveY] = min(derFY(YiZero));
            [valEnterY, idxEnterY] = min(-derFY(YiOne));
            
            if noOut == delta && valLeaveY+valEnterY<valEnterY
                valReduceY = valLeaveY+valEnterY;
                changeModeY = 0;
            elseif noOut == delta
                valReduceY = valEnterY;
                changeModeY = 2;
            elseif valLeaveY<valEnterY
                valReduceY = valLeaveY;
                changeModeY = 1;
            else
                valReduceY = valEnterY;
                changeModeY = 2;
            end
            
            if valReduceY<0
                convergeY = 0;
            else
                convergeY = 1;
            end
        elseif convergeY == 0
            %disp('entering y');
            %             if changeModeY == 0
            %                 %disp('situ 1');
            %                 ri(YiZero(idxLeaveY)) = 1;
            %                 ri(YiOne(idxEnterY)) = 0;
            %                 tmp = YiZero(idxLeaveY);
            %                 YiZero(idxLeaveY) = YiOne(idxEnterY);
            %                 YiOne(idxEnterY) = tmp;
            %
            %             elseif changeModeY == 1
            %                 %disp('situ 2');
            %                 ri(YiZero(idxLeaveY)) = 1;
            %                 tmp = YiZero(idxLeaveY);
            %                 YiZero(idxLeaveY) = [];
            %                 YiOne(end+1) = tmp;
            %                 noOut = noOut+1;
            %             else
            %                 %disp('situ 3');
            %                 ri(YiOne(idxEnterY)) = 0;
            %                 tmp = YiOne(idxEnterY);
            %                 YiOne(idxEnterY) = [];
            %                 YiZero(end+1) = tmp;
            %                 noOut = noOut-1;
            %             end
            % LP2: Fix s and v, solve u
            [ri, ~] = solve_u_lpPR(derFY,delta);          % u can be solved in close form
            derFX(n+1:end) = N'*ri;
            convergeX = 0;
            convergeY = 1;
        end
        
        currF = evalF(ri,si,yi);
        currFVec = [currFVec currF];
        runtimeVec = [runtimeVec toc];
    end
end

%compute final result
theta = computeTheta(yi);
xn = [ri; si];
yn = yi;
currF = evalF(ri,si,yi)

    function [F] = evalF(r, s, y)
        F = abs(r'*(N*y+b) + s'*ones(n,1));
    end

    function tt = computeTheta(yi)
        % From v, compute theta  by substracting the (d+1)-th element from
        % the first d elements of v
        tt = yi(1:d) - repmat(yi(d+1),d,1);
    end
end
