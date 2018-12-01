function [u flag] = solve_u_lpPRNL_Homo(f,delta)
sizeI = length(f)/4;
u = ones(length(f),1);
idxPos = f>=0;
u(idxPos) = 0;
indexNeg = ~idxPos;
idxNeg = find(indexNeg == 1);
if(~isempty(idxNeg))
    fi = zeros(sizeI,1);
    fi(floor(idxNeg/4)+1) = fi(floor(idxNeg/4)+1)+f(idxNeg);
    [fiSorted, idxSorted] = sort(fi,'descend');
    idxSetZero = 1:(sizeI-delta);
    u(idxSorted(idxSetZero)*4) = 0;
    u(idxSorted(idxSetZero)*4-1) = 0;
    u(idxSorted(idxSetZero)*4-2) = 0;
    u(idxSorted(idxSetZero)*4-3) = 0;
end
flag = 1;
end
