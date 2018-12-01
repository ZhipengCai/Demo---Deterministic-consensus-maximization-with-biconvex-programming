%% Generate random data
%  N: number of data points
%  d: dimensions
    
%%
function [x,y,m,c,inlNoise, truthFit, inl,outl, th] = genRandomLinearData(N, dim, in_var, th_in, out_var, outlierP, balance)
 
    if (nargin < 7)
        balance = 1;
    end
    rng('shuffle');
   
    sig = in_var; % Inlier Varience
    osig = out_var;   % Outlier Varience
    th = th_in;  % Inlier Threshold    
    n = dim; % Dimension of space
    %m = randn(n-1, 1);
    m = rand(n-1, 1);
    
    if (dim==2) m = 0.05; disp('hard code m');end
    c = rand; 
    
    %% Generate data
    %xo = randn(n-1,N); 
    xo = -1 + 2*rand(n-1,N);  
    yo = m'*xo + repmat(c,1,N); 

    %% Corrupt data by Gaussian noise
    t = outlierP; 
    t = round(N*t/100); 
    
    x = xo;
    gNoise = -th_in+2*th_in*rand(1,N);
    inlNoise = gNoise;    
    y = yo + gNoise; 
    
     % Add outliers.
    outl = 1:t; 
    inl = 1:N;
    inl(outl) = [];
    th = max(abs(gNoise(inl)));
    
    
    %compute actual number of inliers
    inlIdx = find(abs(y-yo)<=th_in);
    outIdx = find(abs(y-yo)>th_in);
    
    %perturb true inliers to the end
    xo = [xo(:,outIdx), xo(:,inlIdx)];
    x = xo;
    yo = [yo(:,outIdx), yo(:,inlIdx)];
    y = [y(:,outIdx), y(:,inlIdx)];
    
    %% Generate outliers   
    for i=1:t        
        r = y(i) - yo(i);
        k1 = 1; k2 = -2*k1;
        if (~balance) k2 = k1; end;
        if (r > 0)
            while(true)
                outDis = th_in+k1*osig*rand;
                if(abs(outDis) > th_in)
                    break;
                end
            end
            y(i) = yo(i)+outDis;
        else 
            while(true)
                outDis = -th_in+k2*osig*rand;
                if(abs(outDis) > th_in)
                    break;
                end
            end
            y(i) = yo(i)+outDis;
        end
    end
    
    %compute actual number of outliers
    outIdx = find(abs(y-yo)>th_in);
    
    x = x'; 
    y = y';
    x = [x, ones(N, 1)]; 
    
    truthFit = myFitTchebycheff(x(inl,:), y(inl));
    %% Re-estimate inliers
    inl = find(abs(x*truthFit-y)<=th);
    outl=1:N;
    outl(inl) = [];  
    
    
    
end




