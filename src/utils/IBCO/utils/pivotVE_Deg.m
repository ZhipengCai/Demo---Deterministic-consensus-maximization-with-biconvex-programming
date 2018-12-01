%pivoting algorithm that can detect degeneracy
function [activeSetLabelOut, activeSetIdxOut, XOut, aInvAOut, bland] = pivotVE_Deg(A,b,XIn,activeSetLabelIn,activeSetIdxIn,aInvAIn,idxMinLambda)
    [m n] = size(aInvAIn);
    
    bland = 0;
    
    activeSetLabelOut = activeSetLabelIn;
    activeSetIdxOut = activeSetIdxIn;
    
    %efficiency can be improved further
    idxUnactive = find(~activeSetLabelOut);
    %find step size alpha and the next activated constraint
    sq = aInvAIn(idxMinLambda,:)';
    denom = A(:,idxUnactive)'*sq;
    idxFeasible = find(denom<-1e-10);
    domi = b(idxUnactive(idxFeasible))-A(:,idxUnactive(idxFeasible))'*XIn;
    denomFeasible = denom(idxFeasible);
    vecAlpha = domi./denomFeasible;
    [alpha, idxNextActive] = min(vecAlpha);
    if (length(find(vecAlpha == alpha))>1 || alpha ==0)
        bland = 1;
    end
    
    %update X
    XOut=XIn+alpha*sq;
    
    idxNextActiveFinal = idxUnactive(idxFeasible(idxNextActive));   
    
    %update label to active set
    activeSetLabelOut(idxNextActiveFinal) = 1;
    activeSetLabelOut(activeSetIdxOut(idxMinLambda)) = 0;
        
    %update index to active set(do not sort the set)
    leaveIdx = activeSetIdxOut(idxMinLambda);
    activeSetIdxOut(idxMinLambda) = idxNextActiveFinal;    

    
    %rank-one update inverse of active set matrix
    u = A(:,idxNextActiveFinal)-A(:,leaveIdx);
    aInvAOut = rankOneUpdateLU_Ac(aInvAIn,u,idxMinLambda);
    if(bland == 1)
        disp('detect degeneracy');
        [activeSetIdxOut, idxSorted] = sort(activeSetIdxOut);
        aInvAOut = aInvAOut(idxSorted,:);
    end
    
end
