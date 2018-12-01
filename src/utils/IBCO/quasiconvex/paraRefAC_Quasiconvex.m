%frank-wolfe algorithm for homography data using parametric reformulation method
function [xn, yn, theta,  currF] = paraRefAC_Quasiconvex(A, b, cd,dd, x0, y0, delta, config, denThresh) %Note: cd, dd are to enforce strictly cheirality contraints
if (nargin < 9)  denThresh = 0; end;
%prepare data
n = size(A,1);
d = size(A,2);
% Build the matrices
N11 = [-A A*ones(d,1)];
% generate cN matrix similar to N11 for A to enforce c_i^T\theta + d_i
% > 0 to prevent possible divide by 0;
cdp = cd';
cN11 = [-cdp cdp*ones(d,1)];
ddp = dd' - denThresh*ones(numel(dd),1);
si = x0(n+1:end);
yi = y0;

disp(['+++Running PR with delta = ' num2str(delta) '......']);
%if delta = 0 solve the lp produced by setting all ui = 0
if delta == 0
    ri = zeros(n,1); %set ui to zero
    Ays = [-N11 -1*eye(n)];
    Ays = [Ays; -eye(n+d+1)];
    bys = [b; zeros(d+1+n,1)];
    
    % Now add the cheirality constraints, only add if user wants denThresh > 0
    if (denThresh > 0)
        Ays = [Ays; [cN11 zeros(size(cN11,1), n)]];
        bys = [bys; ddp];
        disp(['Adding Constraints CiTx+d > ' num2str(denThresh)] );
    end
    
    fys = [zeros(d+1,1); 1-ri];
    %LP1: fix u, solve s and v
    [ys, ~,exitFlag] = feval(config.lpsolver, fys, Ays, bys);
    
    %syy = [ys1 ys]
    if (exitFlag<1)
        disp(exitFlag);
        disp('LP1 failed');
        xn = [ri; si];
        yn = yi;
        theta = computeTheta(yi);
        currF = evalF(ri, si, yi)
        return;
    end
    %extract s and v from LP1
    yi = ys(1:d+1);
    si = ys(d+2:end);
else
    %setup starting point of u
    ri = x0(1:n);
    maxIter = 1e9;     % Safeguard to prevent infinite loop.
    accuracy = 1e-9;
    m=d+1;
    %%initialize a vertex solution
    %1. fixing s and v, solve u
    fr = N11*yi+b;
    [ri, ~] = solve_u_lpPR_Homo(fr,delta);
    
    %2. fixing u, solve s and v
    Ays = [-N11 -1*eye(n)];
    Ays = [Ays; -eye(n+d+1)];
    bys = [b; zeros(d+1+n,1)];
    % Now add the cheirality constraints, only add if user wants denThresh > 0
    if (denThresh > 0)
        Ays = [Ays; [cN11 zeros(size(cN11,1), n)]];
        bys = [bys; ddp];
        disp(['Adding Constraints CiTx+d > ' num2str(denThresh)] );
    end
    
    fys = [zeros(d+1,1); 1-ri];
    [ys, ~,exitFlag] = feval(config.lpsolver, fys, Ays, bys);
    
    %syy = [ys1 ys]
    if (exitFlag<1)
        disp(exitFlag);
        disp('LP1 failed');
        xn = [ri; si];
        yn = yi;
        theta = computeTheta(yi);
        currF = evalF(ri, si, yi)
        return;
    end
    %extract s and v
    yi = ys(1:d+1);
    si = ys(d+2:end);
    
    %%start active set algorithm
    currF = evalF(ri, si, yi)
    if(currF<1e-9)
        %compute final result
        theta = computeTheta(yi);
        xn = [ri; si];
        yn = yi;
        return;
    end
    %initialize active set and other parameters for set U & (S,V)
    %current constraints
    AX = sparse([eye(n), N11; eye(m+n)]');
    bX = sparse([-b; zeros(m+n,1)]);
    
    Xi = [si;yi];
    YiOne = find(ri==1);
    YiZero = find(ri==0);
    
    %current derivative
    derFX = [ones(n,1); N11'*ri];
    derFY = b+N11*yi;
    
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
    aAX = AX(:,aIndexX);
    if(bland == 1)
        [aAX,aSubIndexX] = lisubmat(aAX,'cols');
        aIndexX = aIndexX(aSubIndexX);
        aIdxX = zeros(length(aIdxX),1);
        aIdxX(aIndexX) = 1;
    end
    %abX = bX(aIdxX);

    %length(si)
    %length(yi)
    %size(aAX)
    %inverse matrix for active set
    aInvAX = inv(aAX);

    %flag showing that last iteration has no improvements
    convergeX = 0;
    %calculate whether there is possible improvements on ui
    if valReduceY<0
        convergeY = 0;
    else
        convergeY = 1;
    end
    
    for iter = 1:maxIter  
        %compute multiplier
        lambdaX = aInvAX*derFX;
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
                [~, idxMinNegLambdaX] = max((lambdaX(idxNegLambdaX).*lambdaX(idxNegLambdaX))./sqNorm);
            else
                idxMinNegLambdaX = 1;
            end
            idxMinLambdaX = idxNegLambdaX(idxMinNegLambdaX);
            leaveLambdaX = lambdaX(idxMinLambdaX);
        end
    
        
        if(convergeX==1 && convergeY==1)
            disp('converged');
            iter
            break;
        end                
        
        if convergeX == 0 && (convergeY==1||leaveLambdaX<valReduceY)
            
            if bland == 0
                %disp('entering x');
                [aIdxX, aIndexX, Xi, aInvAX, bland] = pivotVE_Deg(AX,bX,Xi,aIdxX,aIndexX,aInvAX,idxMinLambdaX);
            else
                %disp('entering bland');
                [aIdxX, aIndexX, Xi, aInvAX, convergeX, bland] = pivotVE_Bland(AX,bX,Xi,aIdxX,aIndexX,aInvAX,derFX,idxNegLambdaX); 
            end
            si = Xi(1:n);
            yi = Xi((n+1):end);
            derFY = b+N11*yi;
            
        elseif convergeY == 0
            %disp('entering y');
            if changeModeY == 0
                %disp('situ 1');
                ri(YiZero(idxLeaveY)) = 1;
                ri(YiOne(idxEnterY)) = 0;
                tmp = YiZero(idxLeaveY);
                YiZero(idxLeaveY) = YiOne(idxEnterY);
                YiOne(idxEnterY) = tmp;
                
            elseif changeModeY == 1
                %disp('situ 2');
                ri(YiZero(idxLeaveY)) = 1;
                tmp = YiZero(idxLeaveY);
                YiZero(idxLeaveY) = [];
                YiOne(end+1) = tmp;
                noOut = noOut+1;
            else
                %disp('situ 3');
                ri(YiOne(idxEnterY)) = 0;
                tmp = YiOne(idxEnterY);
                YiOne(idxEnterY) = [];
                YiZero(end+1) = tmp;
                noOut = noOut-1;
            end
                
            derFX(n+1:end) = N11'*ri;
            convergeX = 0;
            

        end 
        
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
        
    end
    
end

xn = [ri; si];
yn = yi;
theta = computeTheta(yi);
currF = evalF(ri, si, yi)




    function [F] = evalF(r, s, y)
        F = abs(r'*(N11*y+b) + s'*ones(n,1));
    end

    function tt = computeTheta(yi)
        tt = yi(1:d) - repmat(yi(d+1),d,1);
    end

end


