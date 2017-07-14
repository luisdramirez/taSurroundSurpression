%%% analyzeSurroundSurpression
function [rawData] = analyzeSurroundSurpression(subject, runNumber)

subject = 'Pilot';
runNumber = 1;

plotData = 'Yes';

expDir = pwd;
dataDir = 'data';
cd(dataDir)

if exist(['vTA_surrSuppressionOneStim_', subject, '.mat'],'file') ~= 0
    load(['vTA_surrSuppressionOneStim_', subject, '.mat']);
else
    error('Data file does not exist.')
end

targetContrasts = theData(runNumber).p.t1Contrasts; 
estimatedContrast = theData(runNumber).data.estimatedContrast; %subject response/perceived contrast
differenceContrast = theData(runNumber).data.differenceContrast; %how far off response is from target
responseTime = theData(runNumber).data.responseTime;

% [stimConfig t1Contrast estimatedContrast differenceContrast targetOrientation]
rawData = [theData(runNumber).p.trialEvents(:,1),...
    theData(runNumber).p.trialEvents(:,2), theData(runNumber).data.estimatedContrast,...
    theData(runNumber).data.differenceContrast theData(runNumber).p.trialEvents(:,3)];

% Raw data separated into configurations
collTrialsIndx = rawData(:,1) == 1; 
orthTrialsIndx = rawData(:,1) == 2;
baseTrialsIndx = rawData(:,1) == 3;

% [stimConfig t1Contrast estimatedContrast differenceContrast targetOrientation]
collTrials = rawData(collTrialsIndx(:),:); 
orthTrials = rawData(orthTrialsIndx(:),:); 
baseTrials = rawData(baseTrialsIndx(:),:); 

collContrastAvg = zeros(1,length(targetContrasts));
orthContrastAvg = zeros(1,length(targetContrasts));
baseContrastAvg = zeros(1,length(targetContrasts));

% Organize averages by contrast
for nContrast = 1:length(targetContrasts)
   collContrastsAvg(nContrast) = mean(collTrials(collTrials(:,2)==targetContrasts(nContrast),3));
   orthContrastsAvg(nContrast) = mean(orthTrials(orthTrials(:,2)==targetContrasts(nContrast),3)); 
   basedContrastsAvg(nContrast) = mean(baseTrials(baseTrials(:,2)==targetContrasts(nContrast),3)); 
end

contrastAvgs = [collContrastsAvg; orthContrastsAvg; basedContrastsAvg];

%% PLOT DATA
if strcmp(plotData, 'Yes')
     
    plot(targetContrasts, contrastAvgs(1,:))
    hold on
    plot(targetContrasts, contrastAvgs(2,:))
    plot(targetContrasts, contrastAvgs(3,:))
    plot(0:0.1:1,0:0.1:1)
    title('contrast vs. perceived contrast')
    legend('coll','ortho','base')
    xlabel('contrasts')
    ylabel('perceived contrast')
    axis square
    ylim([0 1])
    
end
%%
cd(expDir)
end


