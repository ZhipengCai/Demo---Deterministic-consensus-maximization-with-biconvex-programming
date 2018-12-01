%bland's rule for pivoting algorithm
function [activeSetLabelOut, activeSetIdxOut, XOut, aInvAOut, noImproveFlag, bland] = pivotVE_Bland(A,b,XIn,activeSetLabelIn,activeSetIdxIn,aInvAIn,derF,idxMinusLambda)
    [m n] = size(aInvAIn);
    %initialize searching list
    alpha = 0;
    noImproveFlag = 0;
    numIter = 0;
    
    XOut = XIn;
    activeSetLabelOut = activeSetLabelIn;
    activeSetIdxOut = activeSetIdxIn;
    aInvAOut = aInvAIn;
    

     
    %efficiency can be improved further
    while(alpha == 0)
        bland = 0; 
        idxUnactive = find(~activeSetLabelOut);
        %find step size alpha and the next activated constraint
        sq = aInvAOut(idxMinusLambda(1),:)';
        denom = A(:,idxUnactive)'*sq;
        idxFeasible = find(denom<-1e-9);
        
        domi = b(idxUnactive(idxFeasible))-A(:,idxUnactive(idxFeasible))'*XIn;
        denomFeasible = denom(idxFeasible);
        vecAlpha = domi./denomFeasible;
        [alpha, idxNextActive] = min(vecAlpha);
        if (length(find(vecAlpha == alpha))>1 || alpha==0)
            bland = 1;
        end
        
        idxNextActiveFinal = idxUnactive(idxFeasible(idxNextActive));
        
        %update X
        XOut=XOut+alpha*sq;
        %update label to active set
        activeSetLabelOut(idxNextActiveFinal) = 1;
        activeSetLabelOut(activeSetIdxOut(idxMinusLambda(1))) = 0;
        %update index to active set
        leaveIdx = activeSetIdxOut(idxMinusLambda(1));
        activeSetIdxOut = find(activeSetLabelOut == 1);
        

        idxPermute = find(activeSetIdxOut == idxNextActiveFinal);
        aInvAOut(idxMinusLambda(1),:) = [];
        aInvAOut(idxPermute:n,:) = [sq'; aInvAOut(idxPermute:(n-1),:)]; 

        
        %rank-one update inverse of active set matrix
        u = A(:,idxNextActiveFinal)-A(:,leaveIdx);         
        aInvAOut = rankOneUpdateLU_Ac(aInvAOut,u,idxPermute);
    
        if alpha == 0
            numIter=numIter+1
            lambda = aInvAOut*derF;
            idxMinusLambda = find(lambda<-1e-9);
            if isempty(idxMinusLambda)
                disp('converged pivot using Bland rule');
                noImproveFlag = 1;
                return;
            end
        end

    end
    

end
