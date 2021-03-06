%%% analyzeSurroundSurpression
function [rawData] = analyzeSurroundSurpression(subject)

subject = 'Pre-Pilot_YW';

plotData = 'Yes';

expDir = pwd;
dataDir = 'data';
cd(dataDir)

if exist(['vTA_surrSuppression_', subject, '.mat'],'file') ~= 0
    load(['vTA_surrSuppression_', subject, '.mat']);
    runNumber = length(theData);
else
    error('Data file does not exist.')
end

stimConfigs = theData(runNumber).p.stimConfigurations;
targetContrasts = theData(runNumber).p.t1Contrasts;
estimatedContrast = theData(runNumber).data.estimatedContrast; %subject response
differenceContrast = theData(runNumber).data.differenceContrast; %how far off from target
responseTime = theData(runNumber).data.responseTime;

% [stimConfig cueValidity t1Contrast t2Contrast estimatedContrast differenceContrast targetOrientation]
rawData = [theData(runNumber).p.trialEvents(:,1), theData(runNumber).p.trialEvents(:,end), theData(runNumber).p.trialEvents(:,2),...
    theData(runNumber).p.trialEvents(:,3), estimatedContrast,...
    differenceContrast theData(runNumber).p.trialEvents(:,4)];

% Raw data separated into configurations [t1cued; t2cued]
collTrialsIndx = [rawData(:,1) == 1 rawData(:,1) == 2]; 
orthTrialsIndx = [rawData(:,1) == 3 rawData(:,1) == 4];
baseTrialsIndx = [rawData(:,1) == 5 rawData(:,1) == 6];

% [stimConfig cueValidity t1Contrast t2Contrast estimatedContrast differenceContrast targetOrientation]
collTrials = [rawData(collTrialsIndx(:,1),:); rawData(collTrialsIndx(:,2),:)]; 
orthTrials = [rawData(orthTrialsIndx(:,1),:); rawData(orthTrialsIndx(:,2),:)];
baseTrials = [rawData(baseTrialsIndx(:,1),:); rawData(baseTrialsIndx(:,2),:)];


% Trials separated by condition and target
collT1Trials = [collTrials(collTrials(:,1)==1,2) collTrials(collTrials(:,1)==1,5:6)]; % [cueValidity estimatedContrast differenceContrast] 
collT2Trials = [collTrials(collTrials(:,1)==2,2) collTrials(collTrials(:,1)==2,5:6)];

orthT1Trials = [orthTrials(orthTrials(:,1)==3,2) orthTrials(orthTrials(:,1)==3,5:6)]; 
orthT2Trials = [orthTrials(orthTrials(:,1)==4,2) orthTrials(orthTrials(:,1)==4,5:6)];

baseT1Trials = [baseTrials(baseTrials(:,1)==5,2) baseTrials(baseTrials(:,1)==5,5:6)]; 
baseT2Trials = [baseTrials(baseTrials(:,1)==6,2) baseTrials(baseTrials(:,1)==6,5:6)];

% Trials separated by validity
collValidTrials = collTrials(collTrials(:,2)==1,:);
orthValidTrials = orthTrials(orthTrials(:,2)==1,:);
baseValidTrials = baseTrials(baseTrials(:,2)==1,:);

collInvalidTrials = collTrials(collTrials(:,2)==2,:);
orthInvalidTrials = orthTrials(orthTrials(:,2)==2,:);
baseInvalidTrials = baseTrials(baseTrials(:,2)==2,:);

% Separate data by contrasts for each condition
collValidCuedContrasts = nan(length(collValidTrials),2); % [targetContrast estimatedContrast]
collValidCuedContrastsAvg = zeros(1,length(targetContrasts));
collValidCuedContrastsSTD = zeros(1,length(targetContrasts));

for nTrial = 1:length(collValidTrials)
    if collValidTrials(nTrial,1) == 1 % grab t1 contrast and report
        collValidCuedContrasts(nTrial,1) = collValidTrials(nTrial,3);
        collValidCuedContrasts(nTrial,2) = collValidTrials(nTrial,5);      
    elseif collValidTrials(nTrial,1) == 2 % grab t2 contrast and report
        collValidCuedContrasts(nTrial,1) = collValidTrials(nTrial,4);
        collValidCuedContrasts(nTrial,2) = collValidTrials(nTrial,5);     
    end 
end

orthValidCuedContrasts = nan(length(orthValidTrials),2);
orthValidCuedContrastsAvg = zeros(1,length(targetContrasts));
orthValidCuedContrastsSTD = zeros(1,length(targetContrasts));

for nTrial = 1:length(orthValidTrials)
    if orthValidTrials(nTrial,1) == 3 % grab t1 contrast and report
        orthValidCuedContrasts(nTrial,1) = orthValidTrials(nTrial,3);
        orthValidCuedContrasts(nTrial,2) = orthValidTrials(nTrial,5);
    elseif orthValidTrials(nTrial,1) == 4 % grab t2 contrast and report
        orthValidCuedContrasts(nTrial,1) = orthValidTrials(nTrial,4);
        orthValidCuedContrasts(nTrial,2) = orthValidTrials(nTrial,5);
    end
end

baseValidCuedContrasts = nan(length(baseValidTrials),2);
baseValidCuedContrastsAvg = zeros(1,length(targetContrasts));
baseValidCuedContrastsSTD = zeros(1,length(targetContrasts));

for nTrial = 1:length(baseValidTrials)
    if baseValidTrials(nTrial,1) == 5 % grab t1 contrast and report
        baseValidCuedContrasts(nTrial,1) = baseValidTrials(nTrial,3);
        baseValidCuedContrasts(nTrial,2) = baseValidTrials(nTrial,5);
    elseif baseValidTrials(nTrial,1) == 6 % grab t2 contrast and report
        baseValidCuedContrasts(nTrial,1) = baseValidTrials(nTrial,4);
        baseValidCuedContrasts(nTrial,2) = baseValidTrials(nTrial,5);
    end
end

collInvalidCuedContrasts = nan(length(collInvalidTrials),2);
collInvalidCuedContrastsAvg = zeros(1,length(targetContrasts));
collInvalidCuedContrastsSTD = zeros(1,length(targetContrasts));

for nTrial = 1:length(collInvalidTrials)
    if collInvalidTrials(nTrial,1) == 1 % grab t1 contrast and report
        collInvalidCuedContrasts(nTrial,1) = collInvalidTrials(nTrial,3);
        collInvalidCuedContrasts(nTrial,2) = collInvalidTrials(nTrial,5);
    elseif collInvalidTrials(nTrial,1) == 2 % grab t2 contrast and report
        collInvalidCuedContrasts(nTrial,1) = collInvalidTrials(nTrial,4);
        collInvalidCuedContrasts(nTrial,2) = collInvalidTrials(nTrial,5);
    end

end

orthInvalidCuedContrasts = nan(length(orthInvalidTrials),2);
orthInvalidCuedContrastsAvg = zeros(1,length(targetContrasts));
orthInvalidCuedContrastsSTD = zeros(1,length(targetContrasts));

for nTrial = 1:length(orthInvalidTrials)
    if orthInvalidTrials(nTrial,1) == 3 % grab t1 contrast and report
        orthInvalidCuedContrasts(nTrial,1) = orthInvalidTrials(nTrial,3);
        orthInvalidCuedContrasts(nTrial,2) = orthInvalidTrials(nTrial,5);
    elseif orthInvalidTrials(nTrial,1) == 4 % grab t2 contrast and report
        orthInvalidCuedContrasts(nTrial,1) = orthInvalidTrials(nTrial,4);
        orthInvalidCuedContrasts(nTrial,2) = orthInvalidTrials(nTrial,5);
    end
end

baseInvalidCuedContrasts = nan(length(baseInvalidTrials),2);
baseInvalidCuedContrastsAvg = zeros(1,length(targetContrasts));
baseInvalidCuedContrastsSTD = zeros(1,length(targetContrasts));

for nTrial = 1:length(baseInvalidTrials)
    if baseInvalidTrials(nTrial,1) == 5 % grab t1 contrast and report
        baseInvalidCuedContrasts(nTrial,1) = baseInvalidTrials(nTrial,3);
        baseInvalidCuedContrasts(nTrial,2) = baseInvalidTrials(nTrial,5);
    elseif baseInvalidTrials(nTrial,1) == 6 % grab t2 contrast and report
        baseInvalidCuedContrasts(nTrial,1) = baseInvalidTrials(nTrial,4);
        baseInvalidCuedContrasts(nTrial,2) = baseInvalidTrials(nTrial,5);
    end
end

% Organize averages by contrast
for nContrast = 1:length(targetContrasts)
   collValidCuedContrastsAvg(nContrast) = mean(collValidCuedContrasts(collValidCuedContrasts(:,1)==targetContrasts(nContrast),2));
   orthValidCuedContrastsAvg(nContrast) = mean(orthValidCuedContrasts(orthValidCuedContrasts(:,1)==targetContrasts(nContrast),2)); 
   baseValidCuedContrastsAvg(nContrast) = mean(baseValidCuedContrasts(baseValidCuedContrasts(:,1)==targetContrasts(nContrast),2)); 
   
   collValidCuedContrastsSTD(nContrast) = std(collValidCuedContrasts(collValidCuedContrasts(:,1)==targetContrasts(nContrast),2));
   orthValidCuedContrastsSTD(nContrast) = std(orthValidCuedContrasts(orthValidCuedContrasts(:,1)==targetContrasts(nContrast),2));
   baseValidCuedContrastsSTD(nContrast) = std(baseValidCuedContrasts(baseValidCuedContrasts(:,1)==targetContrasts(nContrast),2));
   
   collInvalidCuedContrastsAvg(nContrast) = mean(collInvalidCuedContrasts(collInvalidCuedContrasts(:,1)==targetContrasts(nContrast),2));
   orthInvalidCuedContrastsAvg(nContrast) = mean(orthInvalidCuedContrasts(orthInvalidCuedContrasts(:,1)==targetContrasts(nContrast),2)); 
   baseInvalidCuedContrastsAvg(nContrast) = mean(baseInvalidCuedContrasts(baseInvalidCuedContrasts(:,1)==targetContrasts(nContrast),2));
   
   collInvalidCuedContrastsSTD(nContrast) = std(collInvalidCuedContrasts(collInvalidCuedContrasts(:,1)==targetContrasts(nContrast),2));
   orthInvalidCuedContrastsSTD(nContrast) = std(orthInvalidCuedContrasts(orthInvalidCuedContrasts(:,1)==targetContrasts(nContrast),2));
   baseInvalidCuedContrastsSTD(nContrast) = std(baseInvalidCuedContrasts(baseInvalidCuedContrasts(:,1)==targetContrasts(nContrast),2));
end

validContrastAvgs = [collValidCuedContrastsAvg; orthValidCuedContrastsAvg; baseValidCuedContrastsAvg];
invalidContrastAvgs = [collInvalidCuedContrastsAvg; orthInvalidCuedContrastsAvg; baseInvalidCuedContrastsAvg];

% standard error
collValidCuedContrastSTE = collValidCuedContrastsSTD/sqrt(length(collValidTrials));
orthValidCuedContrastSTE = orthValidCuedContrastsSTD/sqrt(length(orthValidTrials));
baseValidCuedContrastSTE = baseValidCuedContrastsSTD/sqrt(length(baseValidTrials));

collInvalidCuedContrastSTE = collInvalidCuedContrastsSTD/sqrt(length(collInvalidTrials));
orthInvalidCuedContrastSTE = orthInvalidCuedContrastsSTD/sqrt(length(orthInvalidTrials));
baseInvalidCuedContrastSTE = baseInvalidCuedContrastsSTD/sqrt(length(baseInvalidTrials));


contrastMatrix = nan(length(targetContrasts),length(targetContrasts),2,length(stimConfigs));

for nConfig = 1:length(stimConfigs)
    for nCue = 1:2   
        for iContrast = 1:length(targetContrasts)     
            for jContrast = 1:length(targetContrasts)
                contrastMatrix(iContrast,jContrast,nCue,nConfig) = mean(rawData(rawData(:,1)==nConfig & rawData(:,2)==nCue & rawData(:,3)==targetContrasts(iContrast) & rawData(:,4)==targetContrasts(jContrast),5));
            end
        end
   end
end
% 
for t1Contrast = 1:length(targetContrasts)
    for t2Contrast = 1:length(targetContrasts)
        contrastMatrix1(t1Contrast,t2Contrast) = mean(rawData(rawData(:,1)==1 & rawData(:,2)==1 & rawData(:,3)==targetContrasts(t1Contrast) & rawData(:,4)==targetContrasts(t2Contrast),5));
    end
end
% 
% for iContrast = 1:length(targetContrasts)
%     contrastMatrix1(iContrast,:) = contrastMatrix1(iContrast,:) / targetContrasts(iContrast); 
% end

%% PLOT DATA
if strcmp(plotData, 'Yes')
     
    errorbar(targetContrasts, validContrastAvgs(1,:),collValidCuedContrastSTE) %collinear valid data w/ error
    hold on
    errorbar(targetContrasts, validContrastAvgs(2,:),orthValidCuedContrastSTE) %orthogonal valid data w/ error
    errorbar(targetContrasts, validContrastAvgs(3,:),baseValidCuedContrastSTE) %baseline valid data w/ error
    plot(0:0.1:1,0:0.1:1)
    title('valid')
    legend('coll','ortho','base','unity')
    xlabel('contrasts')
    ylabel('perceived contrast')
    axis square
    ylim([0 1])
    
    
    
    figure
    errorbar(targetContrasts, invalidContrastAvgs(1,:),collInvalidCuedContrastSTE)  %collinear invalid data w/ error
    hold on
    errorbar(targetContrasts, invalidContrastAvgs(2,:),orthInvalidCuedContrastSTE) %orthogonal invalid data w/ error
    errorbar(targetContrasts, invalidContrastAvgs(3,:),baseInvalidCuedContrastSTE) %baseline invalid data w/ error
    plot(0:0.1:1,0:0.1:1)
    title('invalid')
    legend('coll','ortho','base','unity')
    xlabel('contrasts')
    ylabel('perceived contrast')
    axis square
    ylim([0 1])
    
    % colinear valid v invalid
    figure
    errorbar(targetContrasts, validContrastAvgs(1,:),collValidCuedContrastSTE) %collinear valid data w/ error
    hold on
    errorbar(targetContrasts, invalidContrastAvgs(1,:),collInvalidCuedContrastSTE)  %collinear invalid data w/ error
    plot(0:0.1:1,0:0.1:1)
    title('collinear valid v invalid')
    legend('valid','invalid','unity')
    xlabel('contrast')
    ylabel('perceived contrast')
    axis square
    ylim([0 1])
    
    % orthogonal valid v invalid
    figure
    errorbar(targetContrasts, validContrastAvgs(2,:),orthValidCuedContrastSTE) %orthogonal valid data w/ error
    hold on
    errorbar(targetContrasts, invalidContrastAvgs(2,:),orthInvalidCuedContrastSTE) %orthogonal invalid data w/ error
    plot(0:0.1:1,0:0.1:1)
    title('orthogonal valid v invalid')
    legend('valid','invalid','unity')
    xlabel('contrast')
    ylabel('perceived contrast')
    axis square
    ylim([0 1])
    
    % no surround valid v invalid   
    figure
    errorbar(targetContrasts, validContrastAvgs(3,:),baseValidCuedContrastSTE) %baseline valid data w/ error
    hold on
    errorbar(targetContrasts, invalidContrastAvgs(3,:),baseInvalidCuedContrastSTE) %baseline invalid data w/ error
    plot(0:0.1:1,0:0.1:1)
    title('no surround valid v invalid')
    legend('valid','invalid','unity')
    xlabel('contrast')
    ylabel('perceived contrast')
    axis square
    ylim([0 1])
    
    
%     % plot contrast matrix (iContrast, jContrast, nCue, nConfig)
%     for nConfig = 1:length(stimConfigs)
%         for nCue = 1:2   
%             figure
%             heatmap(contrastMatrix(:,:,nCue,nConfig))
%         end
%     end

end

cd(expDir)
end


