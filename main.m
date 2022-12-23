clear; close; clc;

%% Parameters
maxTrials = 25;     % Maximum # of trials
maxIter = 10^3;     % Maximum # of iterations (Stopping criteria)
%maxnFeval = 10^6;   % Maximum # of Cost function evalutaion (Stopping criteria)
eps = realmin;      % Minimum value of Cost function (Stopping criteria)
nBenchmarkFunc = 17; % # of Benchmark Test Functions, it can be 18 at most

empty_Sol.Problem.funName = '';
empty_Sol.Problem.nVar = 10;
empty_Sol.Problem.VarMin = -10;
empty_Sol.Problem.VarMax = 10;
empty_Sol.xmin = [];
empty_Sol.fmin = Inf;
empty_Sol.nFeval = 0;
empty_Sol.nRef = 0;
empty_Sol.nExp = 0;
empty_Sol.nIC = 0;
empty_Sol.nOC = 0;
empty_Sol.nShrink = 0;
empty_Sol.nIter = 0;
empty_Sol.BestCost = [];


%% Main loop

for nVar=100%[2, 3, 10, 30, 60, 90, 120, 180, 240, 360, 480, 1000]
    fileName = ['Result\Results_' num2str(nVar) '.txt'];
    fid = fopen(fileName,'w+');
    maxnFeval = 10000*nVar;   % Maximum # of Cost function evalutaion (Stopping criteria)

    for funId=1:nBenchmarkFunc
        [funName, ~, VarMin, VarMax] = Cost(zeros(1,nVar),funId);
        fprintf('Cost Function : %s, Range : [%4.3f, %4.3f], Dimension : %d\n\n',funName,VarMin,VarMax,nVar);
        fprintf(fid,'Cost Function : %s, Range : [%4.3f, %4.3f], Dimension : %d\n\n',funName,VarMin,VarMax,nVar);

        %% Nelder-Mead
        NMS_Sol=repmat(empty_Sol,maxTrials,1);
        for k=1:maxTrials
            NMS_Sol(k).Problem.funName = funName;
            NMS_Sol(k).Problem.nVar = nVar;
            NMS_Sol(k).Problem.VarMin = VarMin;
            NMS_Sol(k).Problem.VarMax = VarMax;
            tic;
            [NMS_Sol(k).xmin,fmin(k),nFeval(k),nRef(k),nExp(k), nIC(k), nOC(k), nShrink(k), nIter(k),NMS_Sol(k).BestCost] = NMS('Cost', funId, nVar, VarMin, VarMax, maxIter, maxnFeval, eps);
            T(k) = toc;
            NMS_Sol(k).fmin = fmin(k);
            NMS_Sol(k).nFeval = nFeval(k);
            NMS_Sol(k).nRef = nRef(k);
            NMS_Sol(k).nExp = nExp(k);
            NMS_Sol(k).nIC = nIC(k);
            NMS_Sol(k).nOC = nOC(k);
            NMS_Sol(k).nShrink = nShrink(k);
            NMS_Sol(k).nIter = nIter(k);
            NMS_Sol(k).elapsedTime = T(k);
        end

        % Display results
        [~,minInd] = min(fmin);
        [~,maxInd] = max(fmin);
        fprintf('\tNelder-Mead after %d iterations : %d\n\n',NMS_Sol(minInd).nIter);
        fprintf(fid,'\tNelder-Mead after %d iterations : %d\n\n',NMS_Sol(minInd).nIter);
        createReport(NMS_Sol,minInd,maxInd,fmin,T,fid);
        fprintf('\n');
        fprintf(fid,'\n');

        %% Adaptive Nelder-Mead
        ANMS_Sol=repmat(empty_Sol,maxTrials,1);
        for k=1:maxTrials
            ANMS_Sol(k).Problem.funName = funName;
            ANMS_Sol(k).Problem.nVar = nVar;
            ANMS_Sol(k).Problem.VarMin = VarMin;
            ANMS_Sol(k).Problem.VarMax = VarMax;
            tic;
            [ANMS_Sol(k).xmin,fmin(k),nFeval(k),nRef(k),nExp(k), nIC(k), nOC(k), nShrink(k), nIter(k),ANMS_Sol(k).BestCost] = ANMS('Cost', funId, nVar, VarMin, VarMax, maxIter, maxnFeval, eps);
            T(k) = toc;
            ANMS_Sol(k).fmin = fmin(k);
            ANMS_Sol(k).nFeval = nFeval(k);
            ANMS_Sol(k).nRef = nExp(k);
            ANMS_Sol(k).nExp = nExp(k);
            ANMS_Sol(k).nIC = nIC(k);
            ANMS_Sol(k).nOC = nOC(k);
            ANMS_Sol(k).nShrink = nShrink(k);
            ANMS_Sol(k).nIter = nIter(k);
            ANMS_Sol(k).elapsedTime = T(k);
        end

        % Display results
        [~,minInd] = min(fmin);
        [~,maxInd] = max(fmin);
        fprintf('\tAdaptive Nelder-Mead after %d iterations : %d\n\n',ANMS_Sol(minInd).nIter);
        fprintf(fid,'\tAdaptive Nelder-Mead after %d iterations : %d\n\n',ANMS_Sol(minInd).nIter);
        createReport(ANMS_Sol,minInd,maxInd,fmin,T,fid);

        %% Improved version of Adaptive Nelder-Mead
        wANMS_Sol=repmat(empty_Sol,maxTrials,1);
        for k=1:maxTrials
            wANMS_Sol(k).Problem.funName = funName;
            wANMS_Sol(k).Problem.nVar = nVar;
            wANMS_Sol(k).Problem.VarMin = VarMin;
            wANMS_Sol(k).Problem.VarMax = VarMax;
            tic
            [wANMS_Sol(k).xmin,fmin(k),nFeval(k),nRef(k),nExp(k), nIC(k), nOC(k), nShrink(k), nIter(k), wANMS_Sol(k).BestCost] = wANMS('Cost', funId, nVar, VarMin, VarMax, maxIter, maxnFeval, eps);
            T(k) = toc;
            wANMS_Sol(k).fmin = fmin(k);
            wANMS_Sol(k).nFeval = nFeval(k);
            wANMS_Sol(k).nRef = nExp(k);
            wANMS_Sol(k).nExp = nExp(k);
            wANMS_Sol(k).nIC = nIC(k);
            wANMS_Sol(k).nOC = nOC(k);
            wANMS_Sol(k).nShrink = nShrink(k);
            wANMS_Sol(k).nIter = nIter(k);
            wANMS_Sol(k).elapsedTime = T(k);
        end

        % Display results
        [~,minInd] = min(fmin);
        [~,maxInd] = max(fmin);
        fprintf('\n\tAccelerated Adaptive Nelder-Mead  with weighted centroids after %d iterations : %d\n\n',wANMS_Sol(minInd).nIter);
        fprintf(fid,'\n\tAccelerated Adaptive Nelder-Mead  with weighted centroids after %d iterations : %d\n\n',wANMS_Sol(minInd).nIter);
        createReport(wANMS_Sol,minInd,maxInd,fmin,T,fid);
        fprintf('\n-----------------------------------------------------------------------\n\n');
        fprintf(fid,'\n-----------------------------------------------------------------------\n\n');
        % Saving Results
        save(['Result\' funName '_' num2str(nVar) '_Results.mat'],'NMS_Sol','ANMS_Sol','wANMS_Sol');

        %% Dolan Moore Performance Analysis
        algos={'NMS', 'ANMS', 'wANMS'};
        % Dolan Moore Performance Analysis based on Running Time
%        figure
%         minT =[NMS_Sol(:).elapsedTime; ANMS_Sol(:).elapsedTime; wANMS_Sol(:).elapsedTime]';
%         dolanMoore(minT, algos, funId, 0);
%         savefig(['figs\runningTime\f' num2str(funId) '_dolanMore_' num2str(nVar) '.fig']);
%         print(gcf,['figs\runningTime\f' num2str(funId)  '_dolanMore_' num2str(nVar)],'-depsc','-r300')
%         print(gcf,['figs\runningTime\f' num2str(funId)  '_dolanMore_' num2str(nVar)],'-dpng','-r300')

        % Dolan Moore Performance Analysis based on best fitness value
        figure
        minF =[NMS_Sol(:).fmin; ANMS_Sol(:).fmin; wANMS_Sol(:).fmin]';
        dolanMoore(minF, algos, funId, 1);
        savefig(['figs\bestsofar\f' num2str(funId)  '_dolanMore_' num2str(nVar) '.fig']);
        print(gcf,['figs\bestsofar\f' num2str(funId)  '_dolanMore_' num2str(nVar)],'-depsc','-r300')
        print(gcf,['figs\bestsofar\f' num2str(funId)  '_dolanMore_' num2str(nVar)],'-dpng','-r300')
    end
    fprintf('\n');
    fprintf(fid,'\n');
    fclose(fid);
    %% Create Tex Files
    createTexFile(nVar,nBenchmarkFunc);
end