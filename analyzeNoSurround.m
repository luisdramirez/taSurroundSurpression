%%% analyzeSurroundSurpression
function [rawData] = analyzeSurroundSurpression(subject)

subject = 'Pre-Pilot-LR';

plotData = 'Yes';

expDir = pwd;
dataDir = 'data';
cd(dataDir)

if exist(['vTA_noSurround_', subject, '.mat'],'file') ~= 0
    load(['vTA_noSurround_', subject, '.mat']);
    runNumber = length(theData);
else
    error('Data file does not exist.')
end

stimConfigNames = theData(runNumber).p.stimConfigurationsNames;
stimConfigs = theData(runNumber).p.stimConfigurations;
targetContrasts = theData(runNumber).p.t1Contrasts;
estimatedContrast = theData(runNumber).data.estimatedContrast; %subject response
differenceContrast = theData(runNumber).data.differenceContrast; %how far off from target
cueNames = theData(runNumber).p.trialCuesNames;

% [stimConfig cueValidity t1Contrast t2Contrast estimatedContrast differenceContrast]
rawData = [theData(runNumber).p.trialEvents(:,1), theData(runNumber).p.trialEvents(:,end), theData(runNumber).p.trialEvents(:,2),...
    theData(runNumber).p.trialEvents(:,3), estimatedContrast,...
    differenceContrast];

% Raw data separated into configurations [t1cued; t2cued]
trialsIndx = [rawData(:,1) == 5 rawData(:,1) == 6];

% [stimConfig cueValidity t1Contrast t2Contrast estimatedContrast differenceContrast targetOrientation]

trials = [rawData(trialsIndx(:,1),:); rawData(trialsIndx(:,2),:)];


% Trials separated by condition and target

t1Trials = [trials(trials(:,1)==5,2) trials(trials(:,1)==5,5:6)]; % [cueValidity estimatedContrast differenceContrast] 
t2Trials = [trials(trials(:,1)==6,2) trials(trials(:,1)==6,5:6)];

% Trials separated by validity

validTrials = trials(trials(:,2)==1,:);

invalidTrials = trials(trials(:,2)==2,:);

% Separate data by contrasts for each condition

validCuedContrasts = nan(length(validTrials),2);
validCuedContrastsAvg = zeros(1,length(targetContrasts));
validCuedContrastsSTD = zeros(1,length(targetContrasts));

for nTrial = 1:length(validTrials)
    if validTrials(nTrial,1) == 5 % grab t1 contrast and report
        validCuedContrasts(nTrial,1) = validTrials(nTrial,3);
        validCuedContrasts(nTrial,2) = validTrials(nTrial,5);
    elseif validTrials(nTrial,1) == 6 % grab t2 contrast and report
        validCuedContrasts(nTrial,1) = validTrials(nTrial,4);
        validCuedContrasts(nTrial,2) = validTrials(nTrial,5);
    end
end

invalidCuedContrasts = nan(length(invalidTrials),2);
invalidCuedContrastsAvg = zeros(1,length(targetContrasts));
invalidCuedContrastsSTD = zeros(1,length(targetContrasts));

for nTrial = 1:length(invalidTrials)
    if invalidTrials(nTrial,1) == 5 % grab t1 contrast and report
        invalidCuedContrasts(nTrial,1) = invalidTrials(nTrial,3);
        invalidCuedContrasts(nTrial,2) = invalidTrials(nTrial,5);
    elseif invalidTrials(nTrial,1) == 6 % grab t2 contrast and report
        invalidCuedContrasts(nTrial,1) = invalidTrials(nTrial,4);
        invalidCuedContrasts(nTrial,2) = invalidTrials(nTrial,5);
    end
end

% Organize averages by contrast
for nContrast = 1:length(targetContrasts)
   validCuedContrastsAvg(nContrast) = mean(validCuedContrasts(validCuedContrasts(:,1)==targetContrasts(nContrast),2)); 
   validCuedContrastsSTD(nContrast) = std(validCuedContrasts(validCuedContrasts(:,1)==targetContrasts(nContrast),2));

   invalidCuedContrastsAvg(nContrast) = mean(invalidCuedContrasts(invalidCuedContrasts(:,1)==targetContrasts(nContrast),2));
   invalidCuedContrastsSTD(nContrast) = std(invalidCuedContrasts(invalidCuedContrasts(:,1)==targetContrasts(nContrast),2));
end

% standard error

validCuedContrastSTE = validCuedContrastsSTD/sqrt(length(validTrials));
invalidCuedContrastSTE = invalidCuedContrastsSTD/sqrt(length(invalidTrials));


contrastMatrix = nan(length(targetContrasts),length(targetContrasts),2,length(stimConfigs));

for nConfig = 1:length(stimConfigs)
    for nCue = 1:2   
        for iContrast = 1:length(targetContrasts)     
            for jContrast = 1:length(targetContrasts)
                contrastMatrix(iContrast,jContrast,nCue,nConfig) = mean(rawData(rawData(:,1)==stimConfigs(nConfig) & rawData(:,2)==nCue & rawData(:,3)==targetContrasts(iContrast) & rawData(:,4)==targetContrasts(jContrast),5));
            end
        end
   end
end

contrastMeans = nan(length(targetContrasts),length(targetContrasts));

for iContrast = 1:length(targetContrasts)
    for jContrast = 1:length(targetContrasts)
        contrastMeans(iContrast,jContrast) = mean([targetContrasts(iContrast) targetContrasts(jContrast)]);    
    end
end

figure
heatmap(targetContrasts, targetContrasts, contrastMeans)
title('contrast averages')

%% PLOT DATA
if strcmp(plotData, 'Yes')
     
    % no surround valid v invalid   
    figure
    errorbar(targetContrasts, validCuedContrastsAvg,validCuedContrastSTE) %baseline valid data w/ error
    hold on
    errorbar(targetContrasts, invalidCuedContrastsAvg,invalidCuedContrastSTE) %baseline invalid data w/ error
    plot(0:0.1:1,0:0.1:1)
    title('no surround valid v invalid')
    legend('valid','invalid','unity')
    xlabel('contrast')
    ylabel('perceived contrast')
    axis square
    ylim([0 1])
    
    
    % plot contrast matrix (iContrast, jContrast, nCue, nConfig)
    for nConfig = 1:length(stimConfigs)
        for nCue = 1:2   
            figure
            heatmap(targetContrasts, targetContrasts, contrastMatrix(:,:,nCue,nConfig))
            title([stimConfigNames(nConfig) cueNames(nCue)])
            
        end
    end

end

cd(expDir)
end


