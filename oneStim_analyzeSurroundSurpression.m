%%% analyzeSurroundSurpression
function analyzeSurroundSurpression(subject)

subject = 'Pre-Pilot_LR';

plotData = 'Yes';

expDir = pwd;
dataDir = 'data';
cd(dataDir)

if exist(['vTA_surrSuppressionOneStim_', subject, '.mat'],'file') ~= 0
    load(['vTA_surrSuppressionOneStim_', subject, '.mat']);
    runNumber = length(theData);
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

collContrastSTD = zeros(1,length(targetContrasts));
orthContrastSTD = zeros(1,length(targetContrasts));
baseContrastSTD = zeros(1,length(targetContrasts));

% Organize averages by contrast
for nContrast = 1:length(targetContrasts)
   collContrastAvg(nContrast) = mean(collTrials(collTrials(:,2)==targetContrasts(nContrast),3));
   orthContrastAvg(nContrast) = mean(orthTrials(orthTrials(:,2)==targetContrasts(nContrast),3)); 
   baseContrastAvg(nContrast) = mean(baseTrials(baseTrials(:,2)==targetContrasts(nContrast),3)); 
   
   collContrastSTD(nContrast) = std(collTrials(collTrials(:,2)==targetContrasts(nContrast),3));
   orthContrastSTD(nContrast) = std(orthTrials(orthTrials(:,2)==targetContrasts(nContrast),3));
   baseContrastSTD(nContrast) = std(baseTrials(baseTrials(:,2)==targetContrasts(nContrast),3));

end

contrastAvgs = [collContrastAvg; orthContrastAvg; baseContrastAvg];
collContrastSTE = collContrastSTD/sqrt(length(collTrials));
orthContrastSTE = orthContrastSTD/sqrt(length(orthTrials));
baseContrastSTE = baseContrastSTD/sqrt(length(baseTrials));

%% PLOT DATA
if strcmp(plotData, 'Yes')
    figure  
%     plot(targetContrasts, contrastAvgs(1,:)) %colinear data

%     plot(targetContrasts, contrastAvgs(2,:)) %orthogonal data
%     plot(targetContrasts, contrastAvgs(3,:)) %baseline data
%     axis square
    ylim([0 1])
    errorbar(targetContrasts, contrastAvgs(1,:), collContrastSTE) %colinear error
    hold on
    errorbar(targetContrasts, contrastAvgs(2,:), orthContrastSTE) %orthogonal error
    errorbar(targetContrasts, contrastAvgs(3,:), baseContrastSTE) %baseline error
    plot(0:0.1:1,0:0.1:1)
    title('contrast vs. perceived contrast')
    xlabel('contrast')
    ylabel('perceived contrast')
    legend('coll','ortho','base','unity')


    
end

cd(expDir)
end


