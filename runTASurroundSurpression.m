
% How does alertness/arousal/temporal attention influence tuned
% normalization cruves?

% If attention enhances orientation-tuned suppression, measures of tuned
% normalization bandwidth may become narrower when attention is directed
% towards a stimulus.

% Code written by IB & YW, modified by LDR

%% PREPARE AND COLLECT INFO

commandwindow
HideCursor;
echo off
clear all
close all
KbName('UnifyKeyNames');
Screen('Preference', 'SkipSyncTests', 0);

p.subject = '101';
<<<<<<< HEAD
p.repetitions = 2; % has to be a multiple of 2 unique repetitions per run
usePowerMate = 'No';
=======
p.numBlocks = 2; % has to be a multiple of 2 unique repetitions per run
usePowerMate = 'Yes';
>>>>>>> 91a165fde8072fa6b68bf121be5ec070765c9b85

switch usePowerMate
    case 'Yes'
        % Check which devicenumber the powermate is assigned to
        powermate = PsychPowerMate('Open');
        if isempty(powermate)
            error('problem with the powermate');
        end

        % Controls the brightness of the powermate color
        PsychPowerMate('SetBrightness', powermate, 200);
    case 'No'
    otherwise
        error('Option not recognized')
end

% Check which devicenumber the keyboard is assigned to
deviceNumber = 0;
[keyBoardIndices, productNames] = GetKeyboardIndices;
% deviceString = 'Corsair Corsair K95W Gaming Keyboard';
% deviceString = 'Apple Inc. Apple Keyboard';
% deviceString = 'Apple Keyboard';
% deviceString = 'CHICONY USB Keyboard';
deviceString = 'Apple Internal Keyboard / Trackpad';
for i = 1:length(productNames)
    if strcmp(productNames{i}, deviceString)
        deviceNumber = keyBoardIndices(i);
        break;
    end
end
if deviceNumber == 0
    error('No device by that name was detected');
end

% Set directories
expDir = pwd; % Set the experimental directory to the current directory 'pwd'
dataDir = 'data'; % Set the path to a directory called 'data'
t.mySeed = sum(100*clock);
rng(t.mySeed); % Make sure to start with a random seed
t.theDate = datestr(now,'yymmdd'); % Collect today's date
t.timeStamp = datestr(now,'HHMM'); % Timestamp

cd(dataDir);
if exist(['vTA_surrSuppression_', p.subject, '.mat'],'file');
    load(['vTA_surrSuppression_', p.subject, '.mat']);
    runNumber = length(theData)+1;
else
    runNumber = 1;
end
cd(expDir);

%% SCREEN PARAMETERS
screens = Screen('Screens'); % look at available screens
p.screenWidthPixels = Screen('Rect', screens(1));
screenWidth = 36; % 29 cm macbook air, 40 cm trinitron crt, 60 cm Qnix screen
viewDistance = 68; % in cm, ideal distance: 1 cm equals 1 visual degree
visAngle = (2*atan2(screenWidth/2, viewDistance))*(180/pi); % Visual angle of the whole screen
p.pixPerDeg = round(p.screenWidthPixels(3)/visAngle); % pixels per degree visual angle
p.grey = 128;

%% SOUND SETUP
InitializePsychSound(1); % 1 for precise timing

% Open audio device for low-latency output
Fs = 44100;
soundAmp = 1;
reqlatencyclass = 2; % level 2 means take full control over the audio device, even if this causes other sound devices to fail or shutdown
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, Fs, 1); % 1= single-channel
t.cueDur = 0.125; % (s)

cueFreqs = [1046.5 440]; % [high C = target 1, lower A = target 2]
for iTone = 1:numel(cueFreqs)
    tone = MakeBeep(cueFreqs(iTone), t.cueDur, Fs);
    cueTones(iTone,:) = applyEnvelope(tone, Fs);
end

%playSound(pahandle, cueTones(1,:)*soundAmp);
%playSound(pahandle, cueTones(2,:)*soundAmp);


%% GRATING PARAMETERS
p.stimConfigurations = [1 2 3 4]; % [5 6] are baseline (no surround)
p.stimConfigurationsNames = {'colinearT1cued' 'colinearT2cued' 'orthogonalT1cued' 'orthogonalT2cued' 'baselineT1cued' 'baseline T2cued'};

% contrast parameters
p.numContrasts = 4; % 4 for piloting
p.minContrast = 0.1;
p.maxContrast = 0.75;
p.t1Contrasts = 10.^linspace(log10(p.minContrast),log10(p.maxContrast),p.numContrasts);
p.t2Contrasts = p.t1Contrasts;
p.surroundContrast = 0.75; 

% size parameters
p.centerSize = round(1 * p.pixPerDeg);
p.surroundSize = p.centerSize * 3;
p.gapSize = round(0.08 * p.pixPerDeg);
p.outerFixation = round(0.05 * p.pixPerDeg);
p.innerFixation = p.outerFixation/1.5;


%% TRIAL EVENTS
% create matrix with all unique trial events based on number of repetitions
% 1 = which perception condition - center-surround: center probed or surround probed; center only; surround only; 
% 2 = center contrast
% 3 = surround contrast
% 4 = orientation gratings

[F1, F2, F3] = BalanceFactors(p.numBlocks, 0, p.stimConfigurations, p.t1Contrasts, p.t2Contrasts);
% fourth column in TrialEvents is order of appearance for surround or center
% order = repmat(1:length(p.stimConfigurations), [p.numContrasts*(p.numContrasts*p.numBlocks), 1]);
% orderBase = repmat(1:length(p.stimConfigurations), [(p.numContrasts*p.numBlocks), 1]);
baselineConditions = repmat(5:6, [p.numContrasts*p.numBlocks, 1]);
contrasts = repmat([p.t1Contrasts' p.t2Contrasts'], [p.numBlocks 1]);
p.trialEvents = [F1, F2, F3];
p.trialEvents = [p.trialEvents;...
    [baselineConditions(:) [[contrasts(:,1) Shuffle(contrasts(:,2))]; ... % zeros(p.numContrasts*p.numBlocks,1)
    [contrasts(:,1) Shuffle(contrasts(:,2))]]  ]]; %

p.numTrials = size(p.trialEvents,1);
p.numTrialsPerBlock = p.numTrials/p.numBlocks;

% every trial should be a random orientation; 
p.targetsOrientation = randsample(1:180, p.numTrials, true); % each target has the same orientation

% surround orientation dependent on colinear/orthogonal condition
p.surroundOrientation = nan(size(p.targetsOrientation));

for n =1:p.numTrials
   if p.trialEvents(n,1) == 1 || p.trialEvents(n,1) == 2 % surround = target if colinear
       p.surroundOrientation(n) = p.targetsOrientation(n);
   elseif p.trialEvents(n,1) == 3 || p.trialEvents(n,1) == 4 % add 90 degrees if orthogonal trial
       p.surroundOrientation(n) = p.targetsOrientation(n) + 90;
   end
end

whichOrientation =  [p.targetsOrientation' p.surroundOrientation'];
p.trialEvents(:, end+1:end+2) = whichOrientation;

% Assign pre,post cues
p.trialCuesNames = {'Valid' 'Invalid'};
cueValidity = 0.75; % cue validity
p.numValidTrials = cueValidity * p.numTrials; p.numInvalidTrials = p.numTrials - p.numValidTrials; % distribution of valid and invalid trials
valids = ones(p.numValidTrials/(p.numTrials/length(p.stimConfigurations)),1);
invalids = 2*ones(p.numInvalidTrials/(p.numTrials/length(p.stimConfigurations)),1);
trialCuesBase = [valids; invalids]; % combine the valid and invalid vectors
trialCues = repmat(trialCuesBase, p.numTrials/length(trialCuesBase),1); 
p.trialEvents(:,end+1) = trialCues;

p.trialEvents % [stimConfiguration, t1Contrast, t2Contrast, targOrientation, surrOrientation, cueValidity]
p.trialEvents = Shuffle(p.trialEvents); 

% Define parameters for the stimulus
freq = 2;
p.freq = p.centerSize/p.pixPerDeg * freq;
p.freqSurround = p.surroundSize/p.pixPerDeg * freq;
p.orientation = 0;
p.phase = randsample(1:180,p.numTrials*4, true);
p.phase = reshape(p.phase, [p.numTrials 4]);
p.probeContrast = randsample(0.1:0.01:0.9, p.numTrials, true);
% p.orientationChecker = [0 90]; % orientations of the checkerboard
% p.phaseChecker = [0 180]; % phases of the checkerboard

%% TIMING PARAMETERS
t.targetDur = 8/60; % nFramesPerTarget/refrate (s)
t.targetSOA = 18/60; %15/60, 16/60, 18/60 (s)
t.retention = 0.8; % (s)
t.feedbackDur = 0.3; % (s)
t.iti = 1; % (s)
t.startTime = 2; % (s)
t.responseTime = []; % (s)
t.cueTargetSOA = 1; % (s)
t.cueLeadTime = 1; %(s)
t.trialDur = t.cueDur + t.targetSOA + t.cueTargetSOA*2 + t.targetDur*2 + t.retention; % duration of the longest trial
t.trialDur = t.trialDur + t.iti;

jit = 0:0.2:1;
trialJit = Shuffle(repmat(jit,1, ceil(p.numTrials/numel(jit))));
t.trialJit = trialJit(1:p.numTrials);

t.runDur = t.trialDur*p.numTrials + sum(t.trialJit);
trialStartTimes = (0:t.trialDur:t.trialDur*p.numTrials-t.trialDur) + cumsum([0 t.trialJit(1:end-1)]);

targetStartTimes = [];
cueStartTimes = [];

% calculate theoretical timing of events
for n = 1:p.numTrials
    targetTimes = [0 t.targetSOA];
    cueTimes = [targetTimes(1)-t.cueTargetSOA targetTimes(2)+t.cueTargetSOA];
    
    targetTimes = targetTimes + trialStartTimes(n) + t.trialJit(n);
    targetStartTimes = [targetStartTimes; targetTimes'];
    
    cueTimes = cueTimes + trialStartTimes(n) + t.trialJit(n);
    cueStartTimes = [cueStartTimes; cueTimes'];
end

t.targetStartTimes = targetStartTimes;
t.cueStartTimes = cueStartTimes;
t.targetEndTimes = targetStartTimes + t.targetDur;
p.numTargets = numel(targetStartTimes);

% t.flickerTime = 0.2; % (s)
% t.flicker = 0.025; % (s)

%% CREATE STIMULI
% make mask to create circle for the center grating
[x,y] = meshgrid((-p.centerSize/2):(p.centerSize/2)-1, (-p.centerSize/2):(p.centerSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
centerGaussian = zeros(p.centerSize); centerGaussian(eccen <= (p.centerSize/2)) = 1;
% Gaussian = conv2(Gaussian, fspecial('gaussian', p.pixPerDeg, p.pixPerDeg), 'same');

% Make transparency mask for aplha blending the two images
centerTransparencyMask = zeros(p.centerSize); centerTransparencyMask(eccen >= ((p.centerSize)/2)) = 255;

% make mask to create circle for the surround grating
[x,y] = meshgrid((-p.surroundSize/2):(p.surroundSize/2)-1, (-p.surroundSize/2):(p.surroundSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
surroundGaussian = zeros(p.surroundSize); surroundGaussian(eccen <= (p.surroundSize/2)) = 1;

surroundTransparencyMask = zeros(p.surroundSize); surroundTransparencyMask(eccen >= ((p.centerSize)/2)) = 255;

% make unique grating for every trial
[Xc,Yc] = meshgrid(0:(p.centerSize-1), 0:(p.centerSize-1));
[Xs,Ys] = meshgrid(0:(p.surroundSize-1), 0:(p.surroundSize-1));

% % Make checkerboard
% checker1 = square( p.freqSurround*2*pi/p.surroundSize * ( Xs.*sin(p.orientationChecker(1)*(pi/180)) + Ys.*cos(p.orientationChecker(1)*(pi/180)) ) - p.phaseChecker(1) );
% checker2 = square( p.freqSurround*2*pi/p.surroundSize * ( Xs.*sin(p.orientationChecker(2)*(pi/180)) + Ys.*cos(p.orientationChecker(2)*(pi/180)) ) - p.phaseChecker(2) );
% fullChecker = (checker1 .* checker2) .* surroundGaussian;
% fullCheckerNeg = fullChecker*-1;
% fullChecker = fullChecker * (p.grey-1) + p.grey;
% fullCheckerNeg = fullCheckerNeg * (p.grey-1) + p.grey;

% Make actual gratings
centerGrating = NaN(p.numTrials*2, p.centerSize, p.centerSize);
surroundGrating = NaN(p.numTrials*2, p.surroundSize, p.surroundSize);
centerTarget = NaN(p.numTrials*2, p.centerSize, p.centerSize);
% surroundTarget = NaN(p.numTrials*2, p.surroundSize, p.surroundSize);

for n = 1:p.numTrials
    center = (sin(p.freq*2*pi/p.centerSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.phase(n,1)));
    centerGrating(n,:,:) = (center .* centerGaussian);
    
    surround = (sin(p.freqSurround*2*pi/p.surroundSize*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.phase(n,2)));
    surroundGrating(n,:,:) = (surround .* surroundGaussian);
    
    targetCenter = (sin(p.freq*2*pi/p.centerSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.phase(n,3)));
    centerTarget(n,:,:) = (targetCenter .* centerGaussian);
    
%     targetSurround = (sin(p.freqSurround*2*pi/p.surroundSize*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.phase(n,4)));
%     surroundTarget(n,:,:) = (targetSurround .* surroundGaussian);
end

%% WINDOW SETUP
[window,rect] = Screen('OpenWindow', screens(2), p.grey);
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
% load('MyGammaTable.mat');
% Screen('LoadNormalizedGammaTable', window, repmat(gammaTable, [1 3]));
HideCursor;
white = 255; green = [0 255 0];

% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Define coordinates where to draw the stimuli
centerX = rect(3)/2; centerY = rect(4)/2;
% Coordinates for location on left and right side of fixation
patch =  [centerX centerY];

Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);


%% START THE EXPERIMENT
% Draw some text to the screen first outside of the experimental loop:

% Experiment setup
space = zeros(1,256); space(KbName('Space')) = 1;
PsychHID('KbQueueCreate',deviceNumber, space);
PsychHID('KbQueueStart', deviceNumber);


welcomeText = ['Welcome!' '\n' '\n'...
    '' '\n' '\n'...
    'On every trial a center and surrounding stimulus are presented at fixation, which both vary in intensities.' '\n' '\n' ...
    'Your task is to actively maintain a precise representation of both intensity levels. ' '\n' '\n' ...
    'On some trials the center and surround will be presented simultaneous, while on other trials' '\n' '\n' ...
    'they will be presented sequentially' '\n' '\n' ...
    'After a 2 second delay period you are asked to dial the probe to match the intensity' '\n' '\n' ...
    'of the probed center or surrounding stimulus you have in memory as closely as possible.' '\n' '\n' ...
    ' ' '\n' '\n' ...
    'Be sure to always maintain steady fixation on the green dot! ' '\n' '\n' '\n' ...
    'Click the dial to continue.' '\n' '\n' ];

DrawFormattedText(window, welcomeText, 'center', 'center', 255);
Screen('Flip', window);

for n=1:size(cueTones,1)
   playSound(pahandle, cueTones(n,:)*soundAmp);
   WaitSecs(1);
end

while 1
    [pmButton, ~] = PsychPowerMate('Get', powermate);
    if pmButton == 1
        break;
    end
end

% Make sure we can press esc to quit the experiment
startKey = zeros(1,256); startKey(KbName({'ESCAPE'})) = 1;
PsychHID('KbQueueCreate', deviceNumber, startKey);


% Preallocate some variables
data.estimatedContrast = NaN(p.numTrials, 1);
data.differenceContrast = NaN(p.numTrials, 1);
data.responseTime = NaN(p.numTrials, 1);

% Make surround and center textures before trial loop
surroundStimulus = nan(p.numTrials,1);
centerStimulus1 = nan(p.numTrials,1);
centerStimulus2 = nan(p.numTrials,1);

for n = 1:p.numTrials
    surroundTexture(:,:,1) = squeeze(surroundGrating(n,:,:)) * (p.surroundContrast* p.grey ) + p.grey;
    surroundTexture(:,:,2) = surroundTransparencyMask;
    surroundStimulus(n) = Screen('MakeTexture', window, surroundTexture);
    
    centerTexture1(:,:,1) = squeeze(centerGrating(n,:,:)) * ( p.trialEvents(n,2) * p.grey ) + p.grey;
%     centerTexture1(:,:,2) = centerTransparencyMask;
    centerStimulus1(n) = Screen('MakeTexture', window, centerTexture1);
    
    centerTexture2(:,:,1) = squeeze(centerGrating(n,:,:)) * ( p.trialEvents(n,3) * p.grey ) + p.grey;
%     centerTexture2(:,:,2) = centerTransparencyMask;
    centerStimulus2(n) = Screen('MakeTexture', window, centerTexture2);
    
%     checkerMaskPos = Screen('MakeTexture', window, fullChecker);
%     checkerMaskNeg = Screen('MakeTexture', window, fullCheckerNeg);
end

centerMask = Screen('MakeTexture', window, centerTransparencyMask);

%------------%
% TRIAL LOOP %
%------------%

nBlock = 1;
trialTimes = zeros(p.numTrials,6); % [start precue t1 t2 postcue end] 
expStart = GetSecs;

for nTrial = 1:p.numTrials
    trialStart = GetSecs;
    trialTimes(nTrial,1) = trialStart;
    
    if nTrial == 1 || nTrial == p.numTrialsPerBlock*(nBlock-1) + 1

        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation]);
        Screen('Flip', window);
        WaitSecs(t.startTime);
    end
    
    % Draw surroundStimulus if not baseline condition
    if p.trialEvents(nTrial,1) ~= 5 || p.trialEvents(nTrial,1) ~= 6  
        Screen('DrawTexture', window, surroundStimulus(nTrial), [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)), p.trialEvents(nTrial,5))
        Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize, p.gapSize)
        %Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize)
    end
    
    % Draw Fixation
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    WaitSecs(t.cueLeadTime);
    
    % Play pre-cue (odd = high tone (T1); even = low tone (T2))
    if mod(p.trialEvents(nTrial,1),2) ~= 0
        playSound(pahandle, cueTones(1,:)*soundAmp);
    else
        playSound(pahandle, cueTones(2,:)*soundAmp);
    end
    
    % Cue-target SOA
    WaitSecs(t.cueTargetSOA);
    
    % Draw centerStimulus1
    Screen('DrawTexture', window, centerMask, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), patch(2)), p.trialEvents(nTrial,4))
    Screen('DrawTexture', window, centerStimulus1(nTrial), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), patch(2)), p.trialEvents(nTrial,4))
    
    % Draw surroundStimulus if not baseline condition
    if p.trialEvents(nTrial,1) ~= 5 || p.trialEvents(nTrial,1) ~= 6  
        Screen('DrawTexture', window, surroundStimulus(nTrial), [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)), p.trialEvents(nTrial,5))
        Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize, p.gapSize)
        %Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize)
    end
    
    % Draw Fixation
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    trialTimes(nTrial,3) = Screen('Flip', window);
    
    % Stim duration
    WaitSecs(t.targetDur);
       
%     for f = 1:t.flickerTime/(t.flicker*2) % display checkerboard
%             
%         Screen('DrawTexture', window, checkerMaskPos, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)));
%         Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
%         Screen('Flip', window);
%         WaitSecs(t.flicker);
% 
%         Screen('DrawTexture', window, checkerMaskNeg, [],CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)));
%         Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
%         Screen('Flip', window);
%         WaitSecs(t.flicker);
% 
%     end
    
    % Draw surroundStimulus if not baseline condition
    if p.trialEvents(nTrial,1) ~= 5 || p.trialEvents(nTrial,1) ~= 6 
        Screen('DrawTexture', window, surroundStimulus(nTrial), [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)), p.trialEvents(nTrial,5))
        Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize, p.gapSize)
        %Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize)    
    end
     
    % Draw Fixation
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    
    % Cue-target SOA interval
    WaitSecs(t.targetSOA);
   
    % Draw centerStimulus2
    Screen('DrawTexture', window, centerMask, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), patch(2)), p.trialEvents(nTrial,4))
    Screen('DrawTexture', window, centerStimulus2(n), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), patch(2)), p.trialEvents(nTrial,4))
    
    % Draw surroundStimulus if not baseline condition
    if p.trialEvents(nTrial,1) ~= 5 || p.trialEvents(nTrial,1) ~= 6  
        Screen('DrawTexture', window, surroundStimulus(nTrial), [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)), p.trialEvents(nTrial,5))
        Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize, p.gapSize)
        %Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize)
    end
    
    % Draw Fixation
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    trialTimes(nTrial,4) = Screen('Flip', window);
    
    % Stim duration
    WaitSecs(t.targetDur);
    
%     for f = 1: t.flickerTime/(t.flicker*2) % display checkerboard
%             
%         Screen('DrawTexture', window, checkerMaskPos, [],CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)));
%         Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
%         Screen('Flip', window);
%         WaitSecs(t.flicker);
% 
%         Screen('DrawTexture', window, checkerMaskNeg, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)));
%         Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
%         Screen('Flip', window);
%         WaitSecs(t.flicker);
% 
%     end
    
    % Draw surroundStimulus if not baseline condition
    if p.trialEvents(nTrial,1) ~= 5 || p.trialEvents(nTrial,1) ~= 6 
        Screen('DrawTexture', window, surroundStimulus(nTrial), [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)), p.trialEvents(nTrial,5))
        Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize, p.gapSize)
        %Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize)
    end
    
    % Draw Fixation
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    
    % cue-target SOA interval
    WaitSecs(t.cueTargetSOA);       
    
    % Play post-cue (odd = high tone (T1); even = low tone (T2))
    if mod(p.trialEvents(nTrial,1),2) ~= 0 && p.trialEvents(nTrial,6) == 1 %if odd and valid
        playSound(pahandle, cueTones(1,:)*soundAmp);
    elseif mod(p.trialEvents(nTrial,1),2) ~= 0 && p.trialEvents(nTrial,6) == 2 %if odd and invalid
        playSound(pahandle, cueTones(2,:)*soundAmp);
    elseif mod(p.trialEvents(nTrial,1),2) == 0 && p.trialEvents(nTrial,6) == 1 %if even and valid
        playSound(pahandle, cueTones(2,:)*soundAmp);
    elseif mod(p.trialEvents(nTrial,1),2) == 0 && p.trialEvents(nTrial,6) == 2 %if even and invalid
        playSound(pahandle, cueTones(1,:)*soundAmp);
    end
    
    % Draw surroundStimulus if not baseline condition
    if p.trialEvents(nTrial,1) ~= 5 || p.trialEvents(nTrial,1) ~= 6 
        Screen('DrawTexture', window, surroundStimulus(nTrial), [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)), p.trialEvents(nTrial,5))
        Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize, p.gapSize)
%         Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize)
    end

    % Draw Fixation
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    trialTimes(nTrial,end) = Screen('Flip', window);
    WaitSecs(t.retention);
    
    % Set up button press
    PsychHID('KbQueueStart', deviceNumber);
    PsychHID('KbQueueFlush');
       
    % Report: show center Grating and allow user to change contrast
    target = Screen('MakeTexture', window, squeeze(centerTarget(nTrial,:,:))* (p.probeContrast(nTrial)*p.grey) + p.grey);
    Screen('DrawTexture', window, target, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), patch(2)), p.trialEvents(nTrial,4))
    
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    startTrial = GetSecs; % get the start time of each trial
    
    switch usePowerMate
        case 'Yes'
            [~, startAngle] = PsychPowerMate('Get', powermate);
            estContrast = p.probeContrast(nTrial);

            while 1 % start inf loop
                % Query PowerMate button state and rotation angle in "clicks"
                [pmButton, angle] = PsychPowerMate('Get', powermate);
                % 1st button is the "or" of the 1st mouse button and the actual PowerMate button
                if startAngle ~= angle

                    % Convert turn of dial first to degrees and then to contrast:
                    %             angles = (angle * 3.8298)/360;
                    angles = ((startAngle-angle)*3.8298);
                    changeContrast = angles/360;
                    estContrast = estContrast - changeContrast; % update the contrast relative to last dial position
                    % Make sure we stay in range

                    if estContrast > 1
                        estContrast = 1;
                    elseif estContrast < 0
                        estContrast = 0.001;
                    end
                    
                    % Show center Grating
                    target = Screen('MakeTexture', window, squeeze(centerTarget(nTrial,:,:))* (estContrast*p.grey) + p.grey);
                    Screen('DrawTexture', window, target, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), patch(2)), p.trialEvents(nTrial,4))

                    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
                    Screen('Flip', window);

                    startAngle = angle;
                end
                if pmButton == 1;
                    data.estimatedContrast(nTrial) = estContrast;
                    if p.trialEvents(nTrial,1) == 1 || p.trialEvents(nTrial,1) == 3
                        data.differenceContrast(nTrial) = p.trialEvents(nTrial,2) - data.estimatedContrast(nTrial);
                    else
                        data.differenceContrast(nTrial) = p.trialEvents(nTrial,3) - data.estimatedContrast(nTrial);
                    end
                    data.responseTime(nTrial) = (GetSecs - startTrial);
                    pmButton = 0;
                    break;
                end
            end
                % check if esc button has been pressed
                [keyIsDown, keyCode] = PsychHID('KbQueueCheck', deviceNumber); %check response
                key = find(keyCode);
                if key == KbName('ESCAPE') % windows = 'esc', mac = 'ESCAPE' If user presses ESCAPE, exit the program.
                    Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
                    Screen('CloseAll');
                    ListenChar(1); % % Go back to unsuppressed mode
                    FlushEvents('keyDown', deviceNumber);
                    error('User exited program.');
                end
            case 'No'
                %Clicks / scroll
                
                % check if esc button has been pressed
                [keyIsDown, keyCode] = PsychHID('KbQueueCheck', deviceNumber); %check response
                key = find(keyCode);
                if key == KbName('ESCAPE') % windows = 'esc', mac = 'ESCAPE' If user presses ESCAPE, exit the program.
                    Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
                    Screen('CloseAll');
                    ListenChar(1); % % Go back to unsuppressed mode
                    FlushEvents('keyDown', deviceNumber);
                    error('User exited program.');
                end
            otherwise
                error('Option not recognized.')
    end

    trialEnd = GetSecs;
    
    trialTimes(nTrial,end) = trialEnd;
    
    % Present center fixation; get ready for next trial
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation]);
    Screen('Flip', window);
    
    WaitSecs(t.iti);
 
    %%% Rest period
    if  nTrial == p.numTrialsPerBlock*nBlock
        rest = GetSecs;
        
        restText = ['Block ' num2str(nBlock) ' of ' num2str(p.numBlocks) ' completed! You can take a short break now, ' '' '\n' ...
            'or press the dial to continue' '\n' '\n' ];
        DrawFormattedText(window, restText, 'center', 'center', white);
        Screen('Flip', window);
        
        nBlock = nBlock+1; 
        
        pmButtonBreak = 0;
        
        while 1
            [pmButtonBreak, a] = PsychPowerMate('Get', powermate);
            if pmButtonBreak == 1;
                break;
            end
        end
        
        t.restTime = (GetSecs-rest)/60;      
    end
 
    
end

t.endTime = (GetSecs-expStart)/60; %Get endtime of the experiment in seconds
%Draw some more text to the screen outside of the loop:
Screen(window,'TextSize',30);
byeByeText = 'Great work! You have finished this run.';
DrawFormattedText(window, byeByeText, 'center', 'center', [255 255 255]);
Screen('Flip', window);
WaitSecs(2);
Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
Screen('CloseAll')

%% SAVE OUT THE DATA FILE
% cd(dataDir);
% theData(runNumber).t = t;
% theData(runNumber).p = p;
% theData(runNumber).data = data;
% eval(['save vTA_surrSuppression_', p.subject, '.mat theData'])
% 
% cd(expDir);
