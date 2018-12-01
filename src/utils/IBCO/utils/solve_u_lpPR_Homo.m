function [u flag] = solve_u_lpPR_Homo(f,delta)
    [fSorted idxSorted] = sort(f,'ascend');
    u = [zeros(length(f),1)];
    idxIJ = idxSorted(fSorted(1:delta)<0); %first delta constraints that are less than 0
    u(idxIJ)=1;
%     %besides setting the first delta uij's to be 1, we can still set all
%     %the uij's with negative residuals and its index i are included in some
%     %of the uij's set to 1 previously
%     idxI = unique(ceil(idxIJ/4)); %corresponding index i in idxIJ
%     %find constraints that are less than 0
%     sizeI = length(idxI);
%     idxTmp = 1:4*sizeI;
%     b = repmat([0; 1; 2; 3],sizeI,1);
%     size(b)
%     size(idxI)
%     idxFN(idxTmp) = idxI(ceil(idxTmp/4))+b;
%     u(idxFN(f(idxFN)<0)) = 1;
     flag = 1;
end