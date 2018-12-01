% Generate constraint matrix Ax<=b from the linear fitting dataset.
% Each data point is associated with two linear constraints 

function [A,b] = genLinearMaxFSMatrix(X, Y, th)
    n = size(X,1);
    %d = size(X,2);
    A=[X;-X];
    b=[Y+repmat(th,n,1); -Y+repmat(th,n,1)];      

    
end