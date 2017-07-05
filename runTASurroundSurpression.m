
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

p.repetitions = 2; % has to be a multiple of 2 unique repetitions per run

% Check which devicenumber the powermate is assigned to
powermate = PsychPowerMate('Open');
if isempty(powermate)
    error('problem with the powermate');
end

% Controls the brightness of the powermate color
PsychPowerMate('SetBrightness', powermate, 20);

% Check which devicenumber the keyboard is assigned to
deviceNumber = 0;
[keyBoardIndices, productNames] = GetKeyboardIndices;
deviceString = 'Corsair Corsair K95W Gaming Keyboard';
% deviceString = 'Apple Inc. Apple Keyboard';
% deviceString = 'Apple Keyboard';
% deviceString = 'CHICONY USB Keyboard';
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
Screens = Screen('Screens'); % look at available screens
p.screenWidthPixels = Screen('Rect', Screens(1));
screenWidth = 36; % 29 cm macbook air, 40 cm trinitron crt, 60 cm Qnix screen
viewDistance = 68; % in cm, ideal distance: 1 cm equals 1 visual degree
visAngle = (2*atan2(screenWidth/2, viewDistance))*(180/pi); % Visual angle of the whole screen
p.pixPerDeg = round(p.screenWidthPixels(3)/visAngle); % pixels per degree visual angle
p.grey = 128;

%% TIMING PARAMETERS
t.stimDur = 1; % (s)
t.retention = 0.8; % Different retention intervalst.iti = 0.3;
t.iti = 0.3; % (s)
t.startTime = 2; % (s)
t.responseTime = []; % (s)
t.flickerTime = 0.2; % (s)
t.flicker = 0.025; % (s)
t.cueDur = 0.1; % (s)

%% SOUND SETUP
InitializePsychSound(1); % 1 for precise timing

% Open audio device for low-latency output
Fs = 44100;
soundAmp = 1;
reqlatencyclass = 2; % level 2 means take full control over the audio device, even if this causes other sound devices to fail or shutdown
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, Fs, 1); % 1= single-channel

cueFreqs = [1046.5 440]; % [high C = target 1, lower A = target 2]
for iTone = 1:numel(cueFreqs)
    tone = MakeBeep(cueFreqs(iTone), t.cueDur, Fs);
    cueTones(iTone,:) = applyEnvelope(tone, Fs);
end

% playSound(pahandle, cueTones(1)*soundAmp)
%% GRATING PARAMETERS
p.stimConfigurations = 4; % 1 = simultaneous, 2 = sequentially - center first, 3 = sequentially - surround first, 4 = ?
p.numContrasts = 1;
p.testContrasts = 10.^linspace(log10(0.1),log10(0.75),p.numContrasts);
p.surroundContrast = p.testContrasts; % same possible surround contrasts as center

%% TRIAL EVENTS
% create matrix with all unique trial events based on number of repetitions
% 1 = which perception condition - center-surround: center probed or surround probed; center only; surround only; 
% 2 = center contrast
% 3 = surround contrast
% 4 = orientation gratings


[F1, F2, F3] = BalanceFactors(p.repetitions, 0, 1:p.stimConfigurations/2, p.testContrasts, p.surroundContrast);
% fourth column in TrialEvents is order of appearance for surround or center
order = repmat(1:2, [p.numContrasts*(p.numContrasts*p.repetitions), 1]);
orderBase = repmat(1:2, [(p.numContrasts*p.repetitions), 1]);
baselineConditions = repmat(3:4, [p.numContrasts*p.repetitions, 1]);
contrasts = repmat([p.testContrasts' p.surroundContrast'], [p.repetitions 1]);
p.trialEvents = [F1, F2, F3];
p.trialEvents = [p.trialEvents order(:); [baselineConditions(:) [[contrasts(:,1) zeros(p.numContrasts*p.repetitions,1)]; ...
    [zeros(p.numContrasts*p.repetitions,1) contrasts(:,2)]] orderBase(:) ]];
% every trial should be a random orientation
p.numTrials = size(p.trialEvents,1); % multiple of locations and possible targets
whichOrientation =  randsample(1:180, p.numTrials, true);
p.trialEvents(:,5) = whichOrientation';
% p.trialEvents = Shuffle(p.trialEvents); % [stimConfiguration, testContrast, surroundContrast, , orientation]

% size parameters
p.centerSize = round(1 * p.pixPerDeg);
p.surroundSize = p.centerSize * 3;
p.gapSize = round(0.08 * p.pixPerDeg);
p.outerFixation = round(0.05*p.pixPerDeg);
p.innerFixation = p.outerFixation/1.5;

% Define parameters for the stimulus
freq = 2;
p.freq = p.centerSize/p.pixPerDeg * freq;
p.freqSurround = p.surroundSize/p.pixPerDeg * freq;
p.orientation = 0;
p.phase = randsample(1:180,p.numTrials*4, true);
p.phase = reshape(p.phase, [p.numTrials 4]);
p.probeContrast = randsample(0.1:0.01:0.9, p.numTrials, true);
p.orientationChecker = [0 90];
p.phaseChecker = [0 180];

%% CREATE STIMULI
% make mask to create circle for the center grating
[x,y] = meshgrid((-p.centerSize/2):(p.centerSize/2)-1, (-p.centerSize/2):(p.centerSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
centerGaussian = zeros(p.centerSize); centerGaussian(eccen <= (p.centerSize/2)) = 1;
% Gaussian = conv2(Gaussian, fspecial('gaussian', p.pixPerDeg, p.pixPerDeg), 'same');

% make mask to create circle for the surround grating
[x,y] = meshgrid((-p.surroundSize/2):(p.surroundSize/2)-1, (-p.surroundSize/2):(p.surroundSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
surroundGaussian = zeros(p.surroundSize); surroundGaussian(eccen <= (p.surroundSize/2)) = 1;

% Make transparency mask for aplha blending the two images
transparencyMask = zeros(p.surroundSize); transparencyMask(eccen >= ((p.centerSize)/2)) = 255;

% make unique grating for every trial
[Xc,Yc] = meshgrid(0:(p.centerSize-1), 0:(p.centerSize-1));
[Xs,Ys] = meshgrid(0:(p.surroundSize-1), 0:(p.surroundSize-1));

% Make checkerboard
checker1 = square( p.freqSurround*2*pi/p.surroundSize * ( Xs.*sin(p.orientationChecker(1)*(pi/180)) + Ys.*cos(p.orientationChecker(1)*(pi/180)) ) - p.phaseChecker(1) );
checker2 = square( p.freqSurround*2*pi/p.surroundSize * ( Xs.*sin(p.orientationChecker(2)*(pi/180)) + Ys.*cos(p.orientationChecker(2)*(pi/180)) ) - p.phaseChecker(2) );
fullChecker = (checker1 .* checker2) .* surroundGaussian;
fullCheckerNeg = fullChecker*-1;
fullChecker = fullChecker * (p.grey-1) + p.grey;
fullCheckerNeg = fullCheckerNeg * (p.grey-1) + p.grey;

% Make actual gratings
centerGrating = NaN(p.numTrials*2, p.centerSize, p.centerSize);
surroundGrating = NaN(p.numTrials*2, p.surroundSize, p.surroundSize);
centerTarget = NaN(p.numTrials*2, p.centerSize, p.centerSize);
surroundTarget = NaN(p.numTrials*2, p.surroundSize, p.surroundSize);

for n = 1:p.numTrials
    center = (sin(p.freq*2*pi/p.centerSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.phase(n,1)));
    centerGrating(n,:,:) = (center .* centerGaussian);
    
    surround = (sin(p.freqSurround*2*pi/p.surroundSize*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.phase(n,2)));
    surroundGrating(n,:,:) = (surround .* surroundGaussian);
    
    targetCenter = (sin(p.freq*2*pi/p.centerSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.phase(n,3)));
    centerTarget(n,:,:) = (targetCenter .* centerGaussian);
    
    targetSurround = (sin(p.freqSurround*2*pi/p.surroundSize*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.phase(n,4)));
    surroundTarget(n,:,:) = (targetSurround .* surroundGaussian);
end

%% WINDOW SETUP
[window,rect] = Screen('OpenWindow', Screens(2), p.grey);
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

while 1
    [pmButton, ~] = PsychPowerMate('Get', powermate);
    if pmButton == 1;
        break;
    end
end

Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation]);
Screen('Flip', window);

startTime = GetSecs;
WaitSecs(t.startTime);

% Make sure we can press esc to quit the experiment
startKey = zeros(1,256); startKey(KbName({'ESCAPE'})) = 1;
PsychHID('KbQueueCreate', deviceNumber, startKey);


% Preallocate some variables
data.estimatedContrast = NaN(p.numTrials, 1);
data.differenceContrast = NaN(p.numTrials, 1);
data.responseTime = NaN(p.numTrials, 1);

%------------%
% TRIAL LOOP %
%------------%

for n = 1:p.numTrials
    
    centerStimulus = Screen('MakeTexture', window, squeeze( centerGrating(n,:,:)) * (p.trialEvents(n,2) * p.grey ) + p.grey);
    
    surroundText(:,:,1) = squeeze( surroundGrating(n,:,:)) * ( p.trialEvents(n,3) * p.grey ) + p.grey;
    surroundText(:,:,2) = transparencyMask;
    surroundStimulus = Screen('MakeTexture', window, surroundText);
    
    checkerMaskPos = Screen('MakeTexture', window, fullChecker);
    checkerMaskNeg = Screen('MakeTexture', window, fullCheckerNeg);
    
    playSound(pahandle, cueTones(1)*soundAmp) % play pre-cue

    if p.trialEvents(n,1) ~= 4 && p.trialEvents(n,4) == 1
        Screen('DrawTexture', window, centerStimulus, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), patch(2)), p.trialEvents(n,5))
    end
    
    if (p.trialEvents(n,1) ~= 3 && p.trialEvents(n,4) == 2) || (p.trialEvents(n,1) == 2 && p.trialEvents(n,4) == 2)
        Screen('DrawTexture', window, surroundStimulus, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)), p.trialEvents(n,5))
        Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize, p.gapSize)
    end
    
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    WaitSecs(t.stimDur);
       
    for f = 1:t.flickerTime/(t.flicker*2)
            
        Screen('DrawTexture', window, checkerMaskPos, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)));
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);

        Screen('DrawTexture', window, checkerMaskNeg, [],CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)));
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);

    end
    
    % Retention interval
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    WaitSecs(t.retention);
    
    % Time to show stimulus in vwm condition    
    if p.trialEvents(n,1) ~= 4 && p.trialEvents(n,4) == 2
        Screen('DrawTexture', window, centerStimulus, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1),patch(2)), p.trialEvents(n,5))
    end
    
    if (p.trialEvents(n,1) ~= 3 && p.trialEvents(n,4) == 1) || (p.trialEvents(n,1) == 2 && p.trialEvents(n,4) == 1)
        Screen('DrawTexture', window, surroundStimulus, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)), p.trialEvents(n,5))
        Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))', p.gapSize, p.gapSize)
    end
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])    
    Screen('Flip', window);
    WaitSecs(t.stimDur);
    
    for f = 1: t.flickerTime/(t.flicker*2)
            
        Screen('DrawTexture', window, checkerMaskPos, [],CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)));
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);

        Screen('DrawTexture', window, checkerMaskNeg, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)));
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);

    end

    playSound(pahandle, cueTones(1)*soundAmp) % play post-cue
    
    % retention interval
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    WaitSecs(t.retention);       
    
    
    % Set up button press
    PsychHID('KbQueueStart', deviceNumber);
    PsychHID('KbQueueFlush');
       
    % Show center or surround Grating and allow user to change contrast
    if p.trialEvents(n,1) == 1 || p.trialEvents(n,1) == 3
        target = Screen('MakeTexture', window, squeeze(centerTarget(n,:,:))* (p.probeContrast(n)*p.grey) + p.grey);
        Screen('DrawTexture', window, target, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), patch(2)), p.trialEvents(n,5))
    elseif p.trialEvents(n,1) == 2 || p.trialEvents(n,1) == 4
        target = Screen('MakeTexture', window, squeeze(surroundTarget(n,:,:))* (p.probeContrast(n)*p.grey) + p.grey);
        Screen('DrawTexture', window, target, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)), p.trialEvents(n,5))
        Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))')
    end
    
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    startTrial = GetSecs; % get the start time of each trial
    
    [~, startAngle] = PsychPowerMate('Get', powermate);
    estContrast = p.probeContrast(n);
    
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
            if p.trialEvents(n,1) == 1 || p.trialEvents(n,1) == 3
                target = Screen('MakeTexture', window, squeeze(centerTarget(n,:,:))* (estContrast*p.grey) + p.grey);
                Screen('DrawTexture', window, target, [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), patch(2)), p.trialEvents(n,5))
            elseif p.trialEvents(n,1) == 2 || p.trialEvents(n,1) == 4
                target = Screen('MakeTexture', window, squeeze(surroundTarget(n,:,:))* (estContrast*p.grey) + p.grey);
                Screen('DrawTexture', window, target, [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], patch(1), patch(2)), p.trialEvents(n,5))
                Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), patch(2))')
            end
            Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
            Screen('Flip', window);
            
            startAngle = angle;
        end
        if pmButton == 1;
            data.estimatedContrast(n) = estContrast;
            if p.trialEvents(n,1) == 1 || p.trialEvents(n,1) == 3
                data.differenceContrast(n) = p.trialEvents(n,2) - data.estimatedContrast(n);
            else
                data.differenceContrast(n) = p.trialEvents(n,3) - data.estimatedContrast(n);
            end
            data.responseTime(n) = (GetSecs - startTrial);
            pmButton = 0;
            break;
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
        
    end
    
    % Present center fixation
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation]);
    Screen('Flip', window);
    WaitSecs(t.iti);
    
    
    %%% Rest period
    if n == ceil(p.numTrials/2)
        rest = GetSecs;
        
        restText = ['Halfway there! You can take a short break now, ' '' '\n' ...
            'or press the dial to continue' '\n' '\n' ];
        DrawFormattedText(window, restText, 'center', 'center', white);
        Screen('Flip', window);
        pmButtonBreak = 0;
        
        while 1
            [pmButtonBreak, a] = PsychPowerMate('Get', powermate);
            if pmButtonBreak == 1;
                break;
            end
        end
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation]);
        Screen('Flip', window);
        WaitSecs(t.iti);
        t.restTime = (GetSecs-rest)/60;
        
    end
    
end

t.endTime = (GetSecs-startTime)/60; %Get endtime of the experiment in seconds
%Draw some more text to the screen outside of the loop:
Screen(window,'TextSize',30);
byeByeText = 'Great work! You have finished this run';
DrawFormattedText(window, byeByeText, 'center', 'center', [255 255 255]);
Screen('Flip', window);
WaitSecs(2);
Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
Screen('CloseAll')

%% SAVE OUT THE DATA FILE
cd(dataDir);
theData(runNumber).t = t;
theData(runNumber).p = p;
theData(runNumber).data = data;
eval(['save vTA_surrSuppression_', p.Subject, '.mat theData'])

cd(expDir);
