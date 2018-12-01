function [pinx, sample_num, T, Tprime] = prosacSample(n, m, qinx, sample_num, psize, T, Tprime)
   
    
    if m == Tprime && sample_num < n
        sample_num = sample_num + 1;
        newT = T*sample_num/(sample_num-psize);
        Tprime = Tprime + ceil(newT - T);
        %T = newT;
    end
   
    % Sample.
    if Tprime >= m
         inx1 = randomsample(qinx(1:(sample_num-1)),(psize-1));
         inx2 = qinx(sample_num);
         pinx = [inx1,inx2];
    else
         pinx = randomsample(qinx(1:sample_num),psize);
    end
    
    %    sample_num_1 = sample_num;
    
               
end