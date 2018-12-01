%selecting enter constriant for pivoting algorithm
function [activeSetLabelOut, activeSetIdxOut, XOut, aInvAOut] = pivotVE(A,b,XIn,activeSetLabelIn,activeSetIdxIn,aInvAIn,idxMinLambda)
    [m n] = size(aInvAIn);

    activeSetLabelOut = activeSetLabelIn;
    activeSetIdxOut = activeSetIdxIn;
    
    %efficiency can be improved further
    idxUnactive = find(~activeSetLabelOut);
    %find step size alpha and the next activated constraint
    sq = aInvAIn(idxMinLambda,:)';
    %AUaT = A(:,idxUnactive)';
    denom = A(:,idxUnactive)'*sq;
    idxFeasible = find(denom<-1e-9);

    domi = b(idxUnactive(idxFeasible))-A(:,idxUnactive(idxFeasible))'*XIn;
    %domi = b(idxUnactive(idxFeasible))-AUaT(idxFeasible,:)*XIn;
    denomFeasible = denom(idxFeasible);
    [alpha, idxNextActive] = min(domi./denomFeasible);

    %update X
    XOut=XIn+alpha*sq;
    
    idxNextActiveFinal = idxUnactive(idxFeasible(idxNextActive));   
    
    %update label to active set
    activeSetLabelOut(idxNextActiveFinal) = 1;
    activeSetLabelOut(activeSetIdxOut(idxMinLambda)) = 0;
%     %update index to active set
%     leaveIdx = activeSetIdxOut(idxMinLambda);
%     activeSetIdxOut = find(activeSetLabelOut == 1);
%     
%     
%     idxPermute = find(activeSetIdxOut == idxNextActiveFinal);
%     aInvAOut(idxMinLambda,:) = [];
%     aInvAOut(idxPermute:n,:) = [sq'; aInvAOut(idxPermute:(n-1),:)];
%     
%     
%     %rank-one update inverse of active set matrix
%     % inv(A+uv') = inv(A)-inv(A)*uv'*inv(A)/1+v'*inv(A)*u
%     u = A(:,idxNextActiveFinal)-A(:,leaveIdx);
%     vT = [zeros(idxPermute-1,1);1;zeros(m-idxPermute,1)]';
%     aInvAOut = rankOneUpdateLU(aInvAOut,u,vT);
        
    %update index to active set(do not sort the set)
    leaveIdx = activeSetIdxOut(idxMinLambda);
    activeSetIdxOut(idxMinLambda) = idxNextActiveFinal;    
    
    %rank-one update inverse of active set matrix
    % inv(A+uv') = inv(A)-inv(A)*uv'*inv(A)/1+v'*inv(A)*u
    % showIdx = [idxNextActiveFinal, leaveIdx]
    u = A(:,idxNextActiveFinal)-A(:,leaveIdx);
    aInvAOut = rankOneUpdateLU_Ac(aInvAIn,u,idxMinLambda);
    
%     %test whether the aInvAOut has computed correctly(is correct!)
%     A
%     aAX = A(:,activeSetIdxOut)
%     aA = inv(aInvAOut)
%     
%     aInvAX = inv(aAX);
%     aInvAX-aInvAOut
%     pause;
end