%%% analyzeSurroundSurpression
function analyzeSurroundSurpression(subject, runNumber)

expDir = pwd;
dataDir = 'data';
cd(dataDir)

if exist(['vTA_surrSuppression_', subject, '.mat'],'file');
    load(['vTA_surrSuppression_', subject, '.mat']);
else
    error('Data file does not exist.')
end

estimatedContrast = theData(runNumber).data.estimatedContrast;
differenceContrast = theData(runNumber).data.differenceContrast;
responseTime = theData(runNumber).data.responseTime;


return