function [u, flag] = solve_u_lpPR(f,delta)
    [fSorted idxSorted] = sort(f,'ascend');
    u = [zeros(length(f),1)];
    u(idxSorted(fSorted(1:delta)<0))=1;
    flag = 1;
end
