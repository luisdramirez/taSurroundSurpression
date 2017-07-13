%%% analyzeSurroundSurpression
function [rawData] = analyzeSurroundSurpression(subject, runNumber)

expDir = pwd;
dataDir = 'data';
cd(dataDir)

if exist(['vTA_surrSuppression_', subject, '.mat'],'file') ~= 0
    load(['vTA_surrSuppression_', subject, '.mat']);
else
    error('Data file does not exist.')
end

estimatedContrast = theData(runNumber).data.estimatedContrast; %subject response
differenceContrast = theData(runNumber).data.differenceContrast; %how far off from target
responseTime = theData(runNumber).data.responseTime;

% [stimConfig cueValidity t1Contrast t2Contrast estimatedContrast differenceContrast targetOrientation]
rawData = [theData(101).p.trialEvents(:,1), theData(101).p.trialEvents(:,end), theData(101).p.trialEvents(:,2), theData(101).p.trialEvents(:,3), ...
    theData(101).data.estimatedContrast, theData(101).data.differenceContrast theData(101).p.trialEvents(:,4)];

% raw data separated into configurations [t1cued; t2cued]
collTrialsIndx = [rawData(:,1) == 1 rawData(:,1) == 2]; 
orthTrialsIndx = [rawData(:,1) == 3 rawData(:,1) == 4];
baseTrialsIndx = [rawData(:,1) == 5 rawData(:,1) == 6];

% [stimConfig cueValidity t1Contrast t2Contrast estimatedContrast differenceContrast targetOrientation]
collTrials = [rawData(collTrialsIndx(:,1),:); rawData(collTrialsIndx(:,2),:)]; 
orthTrials = [rawData(orthTrialsIndx(:,1),:); rawData(orthTrialsIndx(:,2),:)];
baseTrials = [rawData(baseTrialsIndx(:,1),:); rawData(baseTrialsIndx(:,2),:)];

collT1Trials = [collTrials(collTrials(:,1)==1,2) collTrials(collTrials(:,1)==1,5:6)]; % [cueValidity estimatedContrast differenceContrast] 
collT2Trials = [collTrials(collTrials(:,1)==2,2) collTrials(collTrials(:,1)==2,5:6)];

orthT1Trials = [orthTrials(orthTrials(:,1)==3,2) orthTrials(orthTrials(:,1)==3,5:6)]; 
orthT2Trials = [orthTrials(orthTrials(:,1)==4,2) orthTrials(orthTrials(:,1)==4,5:6)];

baseT1Trials = [baseTrials(baseTrials(:,1)==5,2) baseTrials(baseTrials(:,1)==5,5:6)]; 
baseT2Trials = [baseTrials(baseTrials(:,1)==6,2) baseTrials(baseTrials(:,1)==6,5:6)];


% [allT1 allT2; validT1 validT2; invalidT1 invalidT2]
avgCollDiff = [ mean(abs(collT1Trials(:,3))) mean(abs(collT2Trials(:,3)));... % [differenceContrastT1 differenceContrastT2]
    mean(abs( collT1Trials(collT1Trials(:,1)==1,3) )) mean(abs( collT2Trials(collT2Trials(:,1)==1,3) ));...
    mean(abs( collT1Trials(collT1Trials(:,1)==2,3) )) mean(abs( collT2Trials(collT2Trials(:,1)==2,3) ))]; 

avgOrthDiff = [ mean(abs(orthT1Trials(:,3))) mean(abs(orthT2Trials(:,3)));...
    mean(abs( orthT1Trials(orthT1Trials(:,1)==1,3) )) mean(abs( orthT2Trials(orthT2Trials(:,1)==1,3) ));...
    mean(abs( orthT1Trials(orthT1Trials(:,1)==2,3) )) mean(abs( orthT2Trials(orthT2Trials(:,1)==2,3) ))]; 

avgBaseDiff = [ mean(abs(baseT1Trials(:,3))) mean(abs(baseT2Trials(:,3)));...
    mean(abs( baseT1Trials(baseT1Trials(:,1)==1,3) )) mean(abs( baseT2Trials(baseT2Trials(:,1)==1,3) ));...
    mean(abs( baseT1Trials(baseT1Trials(:,1)==2,3) )) mean(abs( baseT2Trials(baseT2Trials(:,1)==2,3) ))];

%[t1Valid t1Invalid; t2Valid t2Invalid]
Y1 = [avgCollDiff(2,1) avgCollDiff(3,1); avgCollDiff(2,2) avgCollDiff(3,2)]; 
bar(Y1)
title('collinear')
legend('valid','invalid')
xlabel('target')
ylabel('average difference contrast')
axis square
ylim([0 1])

Y2 = [avgOrthDiff(2,1) avgOrthDiff(3,1); avgOrthDiff(2,2) avgOrthDiff(3,2)];
figure
bar(Y2)
title('orthogonal')
legend('valid','invalid')
xlabel('target')
ylabel('average difference contrast')
axis square
ylim([0 1])


Y3 = [avgBaseDiff(2,1) avgBaseDiff(3,1); avgBaseDiff(2,2) avgBaseDiff(3,2)];
figure
bar(Y3)
title('baseline')
legend('valid','invalid')
xlabel('target')
ylabel('average difference contrast')
axis square
ylim([0 1])

cd(expDir)
end
