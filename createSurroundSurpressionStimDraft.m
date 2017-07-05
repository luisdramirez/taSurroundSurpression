% createSurroundSurpressionStim Draft

param.orientations = [0 90]; % [0 45 90 135] possible orientations 

param.numTargets = 2; % each target is composed of 2 elements (surround-center)
param.numSurround = param.numTargets; param.numCenter = param.numTargets; % number of center-surround == number of targets
param.numTrials = param.numCenter * param.numTargets * param.numSurround * numel(param.orientations) ;

stimulus.trialOrientations = repmat(param.orientations, ... % [surroundT1 surroundT2 centerT1 centerT2]
    1,...
    2*param.numTargets);

stimulus.trialOrientations = unique(nchoosek(stimulus.trialOrientations,4), 'rows'); % all possible combinations of orientations
stimulus.trialOrientations = Shuffle(stimulus.trialOrientations);

% Grating images for stimuli will be stored here
surroundImgs = cell(size(param.orientations));
centerImgs = cell(size(param.orientations));

% Screen setup
PsychDefaultSetup(2); % default values for PTB
screenNumber = max(Screen('Screens')); % fix display to most external display
param.screenWidthPixels = Screen('Rect', screenNumber); % screen width in pixels 
screenWidth = 36; % cm
viewDistance = 68; %cm
visAngle = (2*atan2(screenWidth/2, viewDistance))*(180/pi); % visual angle of the whole screen
param.pixPerDeg = round(param.screenWidthPixels(3)/visAngle); % pixels per degree visual angle

% Defining colors
param.white = WhiteIndex(screenNumber);
param.grey = white / 2;
param.green = [0 255 0];

% size parameters
param.centerSize = round(1 * param.pixPerDeg);
param.surroundSize = param.centerSize * 3;
param.gapSize = round(0.08 * param.pixPerDeg); % space between inner and outer stimuli
param.outerFixation = round(0.05 * param.pixPerDeg);
param.innerFixation = param.outerFixation / 1.5;

% Parameters for stimulus
freq = 2;
param.freq = param.centerSize / param.pixPerDeg * freq;
param.freqSurround = param.surroundSize / param.pixPerDeg * freq;
param.orientation = 0;
param.phase = randsample(1:180, param.numTrials*4, true);
param.phase = reshape(param.phase, [param.numTrials 4]);
param.probecontrast = randsample(0.1:0.01:0.9, param.numTrials, true);
param.orientationChecker = [0 90];
param.phaseChecker = [0 180];