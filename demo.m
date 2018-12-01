%***********************************************************************
%                      Demo program for paper:
% Z. Cai, T.J. Chin, H. Le and D. Suter
% Deterministic Consensus Maximization with Biconvex Programming, 
% In Proceedings of the European Conference on Computer Vision (ECCV), September, 2018.
%***********************************************************************

%-----------------------------------------------------------------------
%                           WARNING:
% The program is free for non-commercial academic use. Any commercial use is strictly 
% prohibited without the authors' consent. Please acknowledge the authors by citing the 
% above paper in any academic publications that have made use of this package or part of it.
%-----------------------------------------------------------------------

%- If you encounter any problems or questions please email to 
% zhipeng.cai@adelaide.edu.au.
%
% - This demo code uses 'sedumi' as the default linear programming solver.
% If Gurobi was installed on your system with a proper license,
% the solver can be changed to 'gurobi'. Please note that you also
% need to specify the correct path to gurobi by setting the variable
% 'gurobiPath' in the file prepareSolver.m

% - sedumi is slower than Gurobi in general. Gurobi was used as the LP solver
% for all experiments reported in the paper. Note also that different LP solvers
% may return slightly different qualitative results.

% - This demo does not run USAC directly since it is implemented in C++. We
% provide the code of USAC separately in case the reviewers are interested, 
% which can be compliled using cmake in Ubuntu. To run USAC, just enter 'src/util/USAC/build/bin'
% folder, and then run './RunSingleTest' in the command window following
% the instruction provided in the Readme file (provided by the author of USAC).

% - This code was tested on a 64 bit Ubuntu Machine with MATLAB R2017a
% - To use this demo, just run function demo()
%-----------------------------------------------------------------------

function demo()
addpath(genpath('src/'));

solver.LP = 'sedumi'; % solver.LP = 'gurobi' will make gurobi as LP solver, which is generally faster
solver.SOCP = 'sedumi'; % solver for socp
disp('Preparing solver......');
runtimeNote = '';
if ~strcmp(solver.LP,'gurobi') runtimeNote = ' (setting LP solver to Gurobi should make EP and SS faster) '; end
config.solver = prepareSolver(solver);

printSeparator('-');

warning off all

disp('press any key to start demo for linear regression on synthetic data'); pause;
demoLinSynth(); % Demo for synthetic data

disp('press any key to start demo for homography estimation'); pause;
%ratio of correspondences to plot (1 represent to plot all correspondenecs)
corrPlotRateHomo = 1;
% number of runs for each method (recommend to run 50 times to reduce the effect of randomness)
NoRunsHomo = 10;
demoHomo(corrPlotRateHomo, NoRunsHomo); % Demo for real data; ('ave' in the plot means average; 'dev' in the result means standard deviation)
% 
disp('press any key to start demo for fundamental matrix estimation'); pause;
%ratio of correspondences to plot (1 represent to plot all correspondenecs)
corrPlotRateFun = 1;
% number of runs for each method (recommend to run 50 times to reduce the effect of randomness)
NoRunsFun = 10; 
demoFun(corrPlotRateFun, NoRunsFun); % Demo for real data; ('ave' in the plot means average; 'dev' in the result means standard deviation)

    function demoLinSynth()
        printSeparator(' ');printSeparator('*');
        disp(' ROBUST LINEAR FITTING ON SYNTHETIC DATA ');
        printSeparator('*'); printSeparator(' ');
        
        %parameter for data
        N = 1000;   % Number of points
        d = 8;     % data dimension
        outlierP = 50; %outlier rate
        sig = 0.3; % max noise level
        osig = 1.5;  % outlier vanriance
        epsilon = sig;   % Inlier threshold
        
        linearConfig = config;
        linearConfig.QThresh = 1e-12; % numerical accuracy
        
        %parameter for EP
        linearConfig.alpha = 0.5;   % Initial alpha for EP
        linearConfig.kappa = 5;     % Increment rate of alpha for EP
        linearConfig.maxAlpha = 1e10; % max alpha for EP
        linearConfig.gammaSS = 0.01; %smoothing parameter for SS
        
        aveInlsArr = [];
        aveTimeArr = [];
        %-----------------------------------------------------------------------
        printSeparator('-');
        disp(['Generating Random Data with ' num2str(N) ' points, ' num2str(d) ' dimensions and approximately ' num2str(outlierP) '% of outliers...' ]);
        [x,y] = genRandomLinearData(N, d, sig, epsilon, osig, outlierP, 1);
        
        disp('Data generated.');
        %-----------------------------------------------------------------------
        printSeparator('-');
        disp('Executing RANSAC..... ');
        [ransacTheta, ransacInliers, ransacRuntime, ransacSampleSets ] = linearFit(x, y, epsilon, 'RANSAC', randn(d,1), linearConfig);
        disp(['Ransac  #Inliers = '  num2str(ransacInliers)]);
        disp(['Ransac Runtime  = ' num2str(ransacRuntime) ' seconds']);
        disp('Ransac finished.');
        aveInlsArr = [aveInlsArr, ransacInliers];
        aveTimeArr = [aveTimeArr, ransacRuntime];
        %-----------------------------------------------------------------------
        printSeparator('-');
        disp('Executing Lo-RANSAC (LRS)..... ');
        [loransacTheta, loransacInliers, loransacRuntime] = linearFit(x, y, epsilon, 'LRS', randn(d,1), linearConfig, ransacSampleSets);
        disp(['Lo-Ransac  #Inliers = '  num2str(loransacInliers)]);
        disp(['Lo-Ransac Runtime  = ' num2str(loransacRuntime) ' seconds']);
        disp('Lo-Ransac finished.');
        aveInlsArr = [aveInlsArr, loransacInliers];
        aveTimeArr = [aveTimeArr, loransacRuntime];
        %-----------------------------------------------------------------------
        printSeparator('-');
        disp('Executing Fixing Lo-RANSAC (FLRS)..... ');
        [floransacTheta, floransacInliers, floransacRuntime] = linearFit(x, y, epsilon, 'FLRS', randn(d,1), linearConfig, ransacSampleSets);
        disp(['Fixing Lo-Ransac  #Inliers = '  num2str(floransacInliers)]);
        disp(['Fixing Lo-Ransac Runtime  = ' num2str(floransacRuntime) ' seconds']);
        disp('Fixing Lo-Ransac finished.');
        aveInlsArr = [aveInlsArr, floransacInliers];
        aveTimeArr = [aveTimeArr, floransacRuntime];
        %-----------------------------------------------------------------------
        printSeparator('-');
        disp(['Executing the Exact Penalty method (EP) with random initialization ']);
        disp(['inital alpha = ' num2str(linearConfig.alpha)]);
        disp(['kappa = ' num2str(linearConfig.kappa)]);
        [epTheta, epInliers, epRuntime] = linearFit(x, y, epsilon, 'EP', randn(d,1), linearConfig);
        disp(['EP #Inliers = '  num2str(epInliers) ]);
        disp(['EP Runtime = '  num2str(epRuntime) ' seconds' runtimeNote ]);
        disp('EP finished.');
        aveInlsArr = [aveInlsArr, epInliers];
        aveTimeArr = [aveTimeArr, epRuntime];
        %-----------------------------------------------------------------------
        printSeparator('-');
        disp(['Executing the Exact Penalty method (EP) with FLRS results ']);
        disp(['inital alpha = ' num2str(linearConfig.alpha)]);
        disp(['kappa = ' num2str(linearConfig.kappa)]);
        [eprsTheta, eprsInliers, eprsRuntime] = linearFit(x, y, epsilon, 'EP', floransacTheta, linearConfig);
        disp(['FLRS + EP #Inliers = '  num2str(eprsInliers) ]);
        disp(['FLRS + EP Runtime  = ' num2str(floransacRuntime+eprsRuntime) ' seconds' runtimeNote]);
        disp('FLRS + EP finished.');
        aveInlsArr = [aveInlsArr, eprsInliers];
        aveTimeArr = [aveTimeArr, floransacRuntime+eprsRuntime];
        %-----------------------------------------------------------------------
        printSeparator('-');
        disp(['Executing the Smooth Surrogate method (SS) with random initialization ']);
        [ssTheta, ssInliers, ssRuntime] = linearFit(x, y, epsilon, 'SS', randn(d,1), linearConfig);
        disp(['SS #Inliers = '  num2str(ssInliers) ]);
        disp(['SS Runtime = '  num2str(ssRuntime) ' seconds' runtimeNote ]);
        disp('SS finished.');
        aveInlsArr = [aveInlsArr, ssInliers];
        aveTimeArr = [aveTimeArr, ssRuntime];
        %-----------------------------------------------------------------------
        printSeparator('-');
        disp(['Executing the Smooth Surrogate method (SS) with FLRS starting point ']);
        [ssrsTheta, ssrsInliers, ssrsRuntime] = linearFit(x, y, epsilon, 'SS', floransacTheta, linearConfig);
        disp(['SS #Inliers = '  num2str(ssrsInliers) ]);
        disp(['SS + RS Runtime  = ' num2str(floransacRuntime+ssrsRuntime) ' seconds' runtimeNote]);
        disp('SS finished.');
        aveInlsArr = [aveInlsArr, ssrsInliers];
        aveTimeArr = [aveTimeArr, floransacRuntime+ssrsRuntime];
        %-----------------------------------------------------------------------
        printSeparator('-');
        disp(['Executing IBCO with random initialization']);
        [ibcoTheta, ibcoInliers, ibcoRuntime] = linearFit(x, y, epsilon, 'IBCO', randn(d,1), linearConfig);
        disp(['IBCO #Inliers = '  num2str(ibcoInliers) ]);
        disp(['IBCO Runtime = '  num2str(ibcoRuntime) ' seconds']);
        disp('IBCO finished.');
        aveInlsArr = [aveInlsArr, ibcoInliers];
        aveTimeArr = [aveTimeArr, ibcoRuntime];
        %-----------------------------------------------------------------------
        printSeparator('-');
        disp(['Executing IBCO with random initialization']);
        [ibcorsTheta, ibcorsInliers, ibcorsRuntime] = linearFit(x, y, epsilon, 'IBCO', floransacTheta, linearConfig);
        disp(['IBCO #Inliers = '  num2str(ibcorsInliers) ]);
        disp(['FLRS + IBCO Runtime = '  num2str(floransacRuntime + ibcorsRuntime) ' seconds']);
        disp('IBCO finished.');
        aveInlsArr = [aveInlsArr, ibcorsInliers];
        aveTimeArr = [aveTimeArr, floransacRuntime + ibcorsRuntime];
        %-----------------------------------------------------------------------
        printSeparator('=');
        disp('******************** ROBUST LINEAR FITTING RESULTS ***********************');
        disp(['N = ' num2str(N)]);
        disp(['Ransac  #Inliers = '  num2str(ransacInliers)]);
        disp(['Ransac Runtime  = ' num2str(ransacRuntime) ' seconds']);
        printSeparator('.');
        
        disp(['Lo-Ransac  #Inliers = '  num2str(loransacInliers)]);
        disp(['Lo-Ransac Runtime  = ' num2str(loransacRuntime) ' seconds']);
        printSeparator('.');
        
        disp(['Fixing Lo-Ransac  #Inliers = '  num2str(floransacInliers)]);
        disp(['Fixing Lo-Ransac Runtime  = ' num2str(floransacRuntime) ' seconds']);
        printSeparator('.');
        
        disp(['EP #Inliers = '  num2str(epInliers) ]);
        disp(['EP Runtime  = ' num2str(epRuntime) ' seconds ' ]);
        printSeparator('.');
        
        disp(['FLRS + EP #Inliers = '  num2str(eprsInliers) ]);
        disp(['FLRS + EP Runtime  = ' num2str(floransacRuntime+eprsRuntime) ' seconds ' ]);
        printSeparator('.');
        
        disp(['SS #Inliers = '  num2str(ssInliers) ]);
        disp(['SS Runtime  = ' num2str(ssRuntime) ' seconds ' ]);
        printSeparator('.');
        
        disp(['FLRS + SS #Inliers = '  num2str(ssrsInliers) ]);
        disp(['FLRS + SS Runtime  = ' num2str(floransacRuntime+ssrsRuntime) ' seconds ' ]);
        printSeparator('.');
        
        disp(['IBCO #Inliers = '  num2str(ibcoInliers) ]);
        disp(['IBCO Runtime  = ' num2str(ibcoRuntime) ' seconds' ]);
        printSeparator('.');
        
        disp(['FLRS + IBCO #Inliers = '  num2str(ibcorsInliers) ]);
        disp(['FLRS + IBCO Runtime  = ' num2str(floransacRuntime+ibcorsRuntime) ' seconds' ]);
        printSeparator('=');
    end

    function demoHomo(plotRate,noRuns)
        close all;
        printSeparator(' ');printSeparator('*');
        disp(' ROBUST HOMOGRAPHY ESTIMATION ');
        printSeparator('*'); printSeparator(' ');
        epsilon = 4;   % Inlier threshold 
        fconfig = config;
        fconfig.alpha = 0.5; %EP initial penalty parameter 
        fconfig.kappa = 5.0;
        fconfig.gammaSS = 0.01;
        fconfig.noRuns = noRuns;
        fconfig.maxAlpha = 1e8;
        fconfig.QThresh = 1e-12; %accuracy
        
        %method names
        allMethods = {'RANSAC', 'PROSAC', 'GMLE',...% 1 2 3
            'RANSAC-LORANSAC', 'RANSAC-FLORANSAC',...% 4 5
            'EP','FLORANSAC-EP',...% 6 7
            'SS', 'FLORANSAC-SS',...% 8 9
            'IBCO', 'FLORANSAC-IBCO'}; %10 11
        
        
        %legend for result plotting
        allLegends = {'RS', 'PS', 'GMS', 'LRS' , 'FLRS','EP', 'FLRS+EP', 'SS', 'FLRS+SS', 'IBCO', 'FLRS+IBCO'};
        
        runMethodIdx = [1:11]; %select the method to run if u don't want to run through all of them (for method with FLRS initialization, must run RS and FLRS first)
        plotMethodIdx = runMethodIdx; %select the method to plot
        
        methodList = cell(length(runMethodIdx), 1);
        for i = 1:length(runMethodIdx)
            methodList{i} = allMethods{runMethodIdx(i)};
        end
        
        
        fconfig.plotLegend = cell(length(plotMethodIdx));
        for i=1:length(plotMethodIdx)
            plotLegend{i} = allLegends{plotMethodIdx(i)};
        end
        
        datasetFolder = 'data/';
        Dataset = 'NYCBd2';
        
        disp(['Loading data set ', Dataset, ' .....']);
        dataFile = [datasetFolder Dataset '.mat'];
        fData = load(dataFile);
        disp(['finish loading.']);
        data=fData.homographyData{1};
        x1 = data.x1;
        x2 = data.x2;
    
        T1 = data.T1;
        T2 = data.T2;
        matches = data.matches;
        matchingScores = matches.scores;
        th = epsilon*T2(1,1);
        d = 9;
        
        inlMat = zeros(fconfig.noRuns, numel(methodList));
        runtimeMat = zeros(fconfig.noRuns, numel(methodList));
        
        randTheta = rand(d,1);
        for nRun = 1:fconfig.noRuns
            ransacSampleSet = [];
            for mt = 1:numel(methodList)
                mPrevRuntime = 0;
                method = methodList{mt};
                printSeparator('-');
                disp([method '  run ' num2str(nRun)]);
                
                sMethod = strsplit(method, '-');
                
                method=cell2mat(sMethod(end));
                sampleSet = [];
                theta0 = randTheta;
                if numel(sMethod)==2 %Need to read result from previous method
                    if (strcmp(method,'LORANSAC') || strcmp(method, 'FLORANSAC'))
                        sampleSet = ransacSampleSet;
                    else
                        theta0 = FLRSTheta;
                        mPrevRuntime = FLRSRuntime;
                    end
                end
                
                [mTheta, mInls, mRuntime, mSampleSet] = imageHomographyFit(x1, x2, th, method, theta0, sampleSet,matchingScores,fconfig);
                nInls = numel(mInls);
                if(strcmp(method,'RANSAC'))
                    ransacSampleSet = mSampleSet;
                elseif(strcmp(method,'FLORANSAC'))
                    FLRSTheta = mTheta;
                    FLRSRuntime = mRuntime;
                end
                
                if(nRun == 1)
                    if(strcmp(methodList{mt},'FLORANSAC-IBCO'))
                        noInlsToPlot = floor(nInls*plotRate);
                        figure;
                        plot_match(matches, [matches.X1; matches.X2], mInls, 0, noInlsToPlot); %display inliers and outliers found by FLRS+IBCO
                        title(['optimized consensus set (consensus = ',num2str(noInlsToPlot), ') from IBCO'],'fontsize', 11);
                    elseif(strcmp(method,'RANSAC'))
                        figure;
                        plot_match(matches, [matches.X1; matches.X2], [], 1, 0); %display all correspondences
                        title(['input correspondences'],'fontsize', 11);
                        noInlsToPlot = floor(nInls*plotRate);
                        figure;
                        plot_match(matches, [matches.X1; matches.X2], mInls, 0, noInlsToPlot); %display inliers and outliers found by RANSAC
                        title(['optimized consensus set (consensus = ',num2str(noInlsToPlot), ') from RANSAC'],'fontsize', 11);
                    end
                end
                
                %compute final runtime
                mRuntime = mRuntime + mPrevRuntime;
                
                inlMat(nRun,mt) = nInls;
                runtimeMat(nRun,mt) = mRuntime;
                
            end
        end
        
        printSeparator('=');
        disp('***************** HOMOGRAPHY ESTIMATION RESULTS ************************');
        %consensus
        aveInlsArr = zeros(1,numel(methodList));
        devInlsArr = zeros(1,numel(methodList));
        %average runtime
        aveTimeArr = zeros(1,numel(methodList));
        devTimeArr = zeros(1,numel(methodList));
        for i=1:numel(methodList)
            aveInlsArr(i) = mean(inlMat(:,i));
            devInlsArr(i) = std(inlMat(:,i));
            
            aveTimeArr(i) = mean(runtimeMat(:,i));
            devTimeArr(i) = std(runtimeMat(:,i));
            
            disp([methodList{i}, ' Inls = ' num2str(aveInlsArr(i)), ' (dev: ', num2str(devInlsArr(i)),') ',...
                'Time = ', num2str(aveTimeArr(i)), ' (dev: ', num2str(devTimeArr(i)), ') ']);
        end
        
        disp(['number of correspondences: ', num2str(data.N)]);
        printSeparator('=');
        printSeparator('=');
        %plot results
        %consensus
        figure;
        inlBar = bar(aveInlsArr);
        set(gca,'xticklabel',plotLegend,'FontSize',14);
        ylabel('ave consensus', 'Interpreter', 'tex','Fontsize', 18);
        title(['average consensus for homography estimation'],'fontsize', 11);
        grid on;
        
        %deviation consensus (uncomment if u want to see the plot)
        %         figure;
        %         inlDevBar = bar(devInlsArr);
        %         set(gca,'xticklabel',plotLegend,'FontSize',14);
        %         ylabel('consensus deviation', 'Interpreter', 'Latex','Fontsize', 24);
        %         grid on;
        
        %runtime
        figure;
        timeBar = bar(aveTimeArr);
        set(gca,'xticklabel',plotLegend,'FontSize',14);
        ylabel('ave runtime (s)', 'Interpreter', 'tex','Fontsize', 18);
        title(['average runtime for homography estimation'],'fontsize', 11);
        grid on;
    end

    function demoFun(plotRate,noRuns)
        close all;
        printSeparator(' ');printSeparator('*');
        disp(' ROBUST FUNDAMENTAL MATRIX ESTIMATION ');
        printSeparator('*'); printSeparator(' ');
        epsilon = 0.006;   % Inlier threshold (0.002 for Dresden, 0.003 for booksh)
        fconfig = config;
        fconfig.alpha = 1.5; %EP initial penalty parameter (1.5 for Dresen, 2.5 for booksh)
        fconfig.kappa = 1.5;
        fconfig.gammaSS = 0.01;
        fconfig.noRuns = noRuns;
        fconfig.maxAlpha = 1e8;
        fconfig.QThresh = 1e-12; %accuracy
        
        %method names
        allMethods = {'RANSAC', 'PROSAC', 'GMLE',...% 1 2 3
            'RANSAC-LORANSAC', 'RANSAC-FLORANSAC',...% 4 5
            'EP','FLORANSAC-EP',...% 6 7
            'SS', 'FLORANSAC-SS',...% 8 9
            'IBCO', 'FLORANSAC-IBCO'}; %10 11
        
        
        %legend for result plotting
        allLegends = {'RS', 'PS', 'GMS', 'LRS' , 'FLRS','EP', 'FLRS+EP', 'SS', 'FLRS+SS', 'IBCO', 'FLRS+IBCO'};
        
        runMethodIdx = [1:11]; %select the method to run if u don't want to run through all of them (for method 4-7, must run RANSAC first)
        plotMethodIdx = runMethodIdx; %select the method to plot
        
        methodList = cell(length(runMethodIdx), 1);
        for i = 1:length(runMethodIdx)
            methodList{i} = allMethods{runMethodIdx(i)};
        end
        
        
        fconfig.plotLegend = cell(length(plotMethodIdx));
        for i=1:length(plotMethodIdx)
            plotLegend{i} = allLegends{plotMethodIdx(i)};
        end
        
        
        datasetFolder = 'data/';
        Dataset = 'booksh';
        
        disp(['Loading data set ', Dataset, ' .....']);
        dataFile = [datasetFolder Dataset '.mat'];
        fData = load(dataFile);
        disp(['finish loading.']);
        
        data = fData.linearData{1};
        x = data.x;
        y = data.y;
        
        matches = data.matches;
        matchingScores = matches.scores;
        
        d = data.dim;
        
        
        inlMat = zeros(fconfig.noRuns, numel(methodList));
        runtimeMat = zeros(fconfig.noRuns, numel(methodList));
        
        randTheta = rand(d,1);
        for nRun = 1:fconfig.noRuns
            ransacSampleSet = [];
            for mt = 1:numel(methodList)
                mPrevRuntime = 0;
                method = methodList{mt};
                printSeparator('-');
                disp([method '  run ' num2str(nRun)]);
                
                sMethod = strsplit(method, '-');
                
                method=cell2mat(sMethod(end));
                sampleSet = [];
                theta0 = randTheta;
                if numel(sMethod)==2 %Need to read result from previous method
                    if (strcmp(method,'LORANSAC') || strcmp(method, 'FLORANSAC'))
                        sampleSet = ransacSampleSet;
                    else
                        theta0 = FLRSTheta;
                        mPrevRuntime = FLRSRuntime;
                    end
                end
                
                [mTheta, mInls, ~, mRuntime, mSampleSet] = linearFitFundamental(x, y, epsilon, method, theta0, sampleSet, matchingScores,fconfig);
                nInls = numel(mInls);
                if(strcmp(method,'RANSAC'))
                    ransacSampleSet = mSampleSet;
                elseif(strcmp(method,'FLORANSAC'))
                    FLRSTheta = mTheta;
                    FLRSRuntime = mRuntime;
                end
                
                if(nRun == 1)
                    if(strcmp(methodList{mt},'FLORANSAC-IBCO'))
                        noInlsToPlot = floor(nInls*plotRate);
                        figure;
                        plot_match(matches, [matches.X1; matches.X2], mInls, 0, noInlsToPlot); %display inliers and outliers found by IBCO
                        title(['optimized consensus set (consensus = ',num2str(noInlsToPlot), ') from IBCO'],'fontsize', 11);
                    elseif(strcmp(method,'RANSAC'))
                        figure;
                        plot_match(matches, [matches.X1; matches.X2], [], 1, 0); %display all correspondences
                        title(['input correspondences'],'fontsize', 11);
                        noInlsToPlot = floor(nInls*plotRate);
                        figure;
                        plot_match(matches, [matches.X1; matches.X2], mInls, 0, noInlsToPlot); %display inliers and outliers found by RANSAC
                        title(['optimized consensus set (consensus = ',num2str(noInlsToPlot), ') from RANSAC'],'fontsize', 11);
                    end
                end
                
                %compute final runtime
                mRuntime = mRuntime + mPrevRuntime;
                
                inlMat(nRun,mt) = nInls;
                runtimeMat(nRun,mt) = mRuntime;
                
            end
        end
        
        printSeparator('=');
        disp('***************** FUNDAMENTAL MATRIX ESTIMATION RESULTS ************************');
        %consensus
        aveInlsArr = zeros(1,numel(methodList));
        devInlsArr = zeros(1,numel(methodList));
        %average runtime
        aveTimeArr = zeros(1,numel(methodList));
        devTimeArr = zeros(1,numel(methodList));
        for i=1:numel(methodList)
            aveInlsArr(i) = mean(inlMat(:,i));
            devInlsArr(i) = std(inlMat(:,i));
            
            aveTimeArr(i) = mean(runtimeMat(:,i));
            devTimeArr(i) = std(runtimeMat(:,i));
            
            disp([methodList{i}, ' Inls = ' num2str(aveInlsArr(i)), ' (dev: ', num2str(devInlsArr(i)),') ',...
                'Time = ', num2str(aveTimeArr(i)), ' (dev: ', num2str(devTimeArr(i)), ') ']);
        end
        
        disp(['number of correspondences: ', num2str(data.N)]);
        printSeparator('=');
        printSeparator('=');
        %plot results
        %consensus
        figure;
        inlBar = bar(aveInlsArr);
        set(gca,'xticklabel',plotLegend,'FontSize',14);
        ylabel('ave consensus', 'Interpreter', 'tex','Fontsize', 18);
        title(['average consensus for fundamental matrix fitting'],'fontsize', 11);
        grid on;
        
        %deviation consensus (uncomment if u want to see the plot)
        %         figure;
        %         inlDevBar = bar(devInlsArr);
        %         set(gca,'xticklabel',plotLegend,'FontSize',14);
        %         ylabel('consensus deviation', 'Interpreter', 'Latex','Fontsize', 24);
        %         grid on;
        
        %runtime
        figure;
        timeBar = bar(aveTimeArr);
        set(gca,'xticklabel',plotLegend,'FontSize',14);
        ylabel('ave runtime (s)', 'Interpreter', 'tex','Fontsize', 18);
        title(['average runtime for fundamental matrix fitting'],'fontsize', 11);
        grid on;
     end

end
