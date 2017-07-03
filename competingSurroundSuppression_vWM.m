

%%% Do memory respresentations of center and surround interact/compete?

% Measures perceived contrast for either center of surround when both are
% actively maintained in visual working memory

% Code written by IB & YW

%%% PREPARE AND COLLECT INFO
commandwindow
echo off
clear all
close all
KbName('UnifyKeyNames');
Screen('Preference', 'SkipSyncTests', 1);

p.Subject = '004';

p.Repetitions = 2; % has to be a multiple of 2 unique repetitions per run

% Check which devicenumber the powermate is assigned to
powermate = PsychPowerMate('Open');
if isempty(powermate)
    error('problem with the powermate');
end
% Controls the brightness of the powermate color
PsychPowerMate('SetBrightness', powermate, 20);

% Check which devicenumber the keyboard is assigned to
deviceNumber = 0;
[keyBoardIndices, ProductNames] = GetKeyboardIndices;
deviceString = 'Corsair Corsair K95W Gaming Keyboard';
% deviceString = 'Apple Inc. Apple Keyboard';
% deviceString = 'Apple Keyboard';
% deviceString = 'CHICONY USB Keyboard';
for i = 1:length(ProductNames)
    if strcmp(ProductNames{i}, deviceString)
        deviceNumber = keyBoardIndices(i);
        break;
    end
end
if deviceNumber == 0
    error('No device by that name was detected');
end

%set directory
expdir = pwd; %Set the experimental directory to the current directory 'pwd'
datadir = 'data'; %Set the path to a directory called 'Data'
t.MySeed = sum(100*clock);
rng(t.MySeed); % make sure we start with a random seed
t.TheDate = datestr(now,'yymmdd'); %Collect todays date
t.TimeStamp = datestr(now,'HHMM'); %Timestamp

cd(datadir);
if exist(['vWM_Comp_surrSuppression_', p.Subject, '.mat'],'file');
    load(['vWM_Comp_surrSuppression_', p.Subject, '.mat']);
    runnumber = length(TheData)+1;
else
    runnumber = 1;
end
cd(expdir);

%%% SCREEN PARAMETERS
Screens = Screen('Screens'); % look at available screens
p.ScreenWidthPixels = Screen('Rect', Screens(1));
ScreenWidth = 36; % 29 cm macbook air, 40 cm trinitron crt, 60 cm Qnix screen
ViewDistance = 68; % in cm, ideal distance: 1 cm equals 1 visual degree
VisAngle = (2*atan2(ScreenWidth/2, ViewDistance))*(180/pi); % Visual angle of the whole screen
p.ppd = round(p.ScreenWidthPixels(3)/VisAngle); % pixels per degree visual angle
p.Grey = 128;

%%% TIMING PARAMETERS
t.stimon = 1;       % in sec
t.retention = 0.8;    % Different retention intervalst.iti = 0.3;
t.iti = 0.3;
t.Starttime = 2;
t.ResponseTime = [];
t.flickertime = 0.2;
t.flicker = 0.025;

%%% GRATING PARAMETERS
p.stimConfigurations = 4; % 1 = simultaneous, 2 = sequentially - center first, 3 = sequentially - surround first
p.numContrasts = 4;
p.testContrasts = 10.^linspace(log10(0.1),log10(0.75),p.numContrasts);
p.surroundContrast = p.testContrasts; % same possible surround contrasts as center

%%% TRIAL EVENTS
% create matrix with all unique trial events based on number of repetitions
% 1 = which perception condition - center-surround: center probed or surround probed; center only; surround only; 
% 2 = center contrast
% 3 = surround contrast
% 4 = orientation gratings


[F1, F2, F3] = BalanceFactors(p.Repetitions, 0, 1:p.stimConfigurations/2, p.testContrasts, p.surroundContrast);
%fourth column in TrialEvents is order of appearance for surround or center
Order = repmat(1:2, [p.numContrasts*(p.numContrasts*p.Repetitions), 1]);
OrderBase = repmat(1:2, [(p.numContrasts*p.Repetitions), 1]);
BaselineConditions = repmat(3:4, [p.numContrasts*p.Repetitions, 1]);
Contrasts = repmat([p.testContrasts' p.surroundContrast'], [p.Repetitions 1]);
p.TrialEvents = [F1, F2, F3];
p.TrialEvents = [p.TrialEvents Order(:); [BaselineConditions(:) [[Contrasts(:,1) zeros(p.numContrasts*p.Repetitions,1)]; ...
    [zeros(p.numContrasts*p.Repetitions,1) Contrasts(:,2)]] OrderBase(:) ]];
% every trial should be a random orientation
p.numTrials = size(p.TrialEvents,1); % multiple of locations and possible targets
whichOrientation =  randsample(1:180, p.numTrials, true);
p.TrialEvents(:,5) = whichOrientation';
p.TrialEvents = Shuffle(p.TrialEvents);

% size parameters
p.CenterSize = round(1 * p.ppd);
p.SurroundSize = p.CenterSize * 3;
p.GapSize = round(0.08 * p.ppd);
p.OuterFixation = round(0.05*p.ppd);
p.InnerFixation = p.OuterFixation/1.5;

% Define parameters for the stimulus
freq = 2;
p.Freq = p.CenterSize/p.ppd * freq;
p.Freq_surround = p.SurroundSize/p.ppd * freq;
p.orientation = 0;
p.phase = randsample(1:180,p.numTrials*4, true);
p.phase = reshape(p.phase, [p.numTrials 4]);
p.probecontrast = randsample(0.1:0.01:0.9, p.numTrials, true);
p.orientationChecker = [0 90];
p.phaseChecker = [0 180];

%%% CREATE STIMULI
% make mask to create circle for the center grating
[x,y] = meshgrid((-p.CenterSize/2):(p.CenterSize/2)-1, (-p.CenterSize/2):(p.CenterSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
centerGaussian = zeros(p.CenterSize); centerGaussian(eccen <= (p.CenterSize/2)) = 1;
% Gaussian = conv2(Gaussian, fspecial('gaussian', p.ppd, p.ppd), 'same');

% make mask to create circle for the surround grating
[x,y] = meshgrid((-p.SurroundSize/2):(p.SurroundSize/2)-1, (-p.SurroundSize/2):(p.SurroundSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
surroundGaussian = zeros(p.SurroundSize); surroundGaussian(eccen <= (p.SurroundSize/2)) = 1;

% Make transparency mask for aplha blending the two images
transparencyMask = zeros(p.SurroundSize); transparencyMask(eccen >= ((p.CenterSize)/2)) = 255;

% make unique grating for every trial
[Xc,Yc] = meshgrid(0:(p.CenterSize-1),0:(p.CenterSize-1));
[Xs,Ys] = meshgrid(0:(p.SurroundSize-1),0:(p.SurroundSize-1));

% Make checkerboard
checker1 = (square(p.Freq_surround*2*pi/p.SurroundSize*(Xs.*sin(p.orientationChecker(1)*(pi/180))+Ys.*cos(p.orientationChecker(1)*(pi/180)))-p.phaseChecker(1)));
checker2 = (square(p.Freq_surround*2*pi/p.SurroundSize*(Xs.*sin(p.orientationChecker(2)*(pi/180))+Ys.*cos(p.orientationChecker(2)*(pi/180)))-p.phaseChecker(2)));
FullChecker = (checker1.*checker2) .* surroundGaussian;
FullCheckerNeg = (FullChecker*-1);
FullChecker = FullChecker * (p.Grey-1) + p.Grey;
FullCheckerNeg = FullCheckerNeg * (p.Grey-1) + p.Grey;

% Make actual gratings
center_grating = NaN(p.numTrials*2,p.CenterSize,p.CenterSize);
surround_grating = NaN(p.numTrials*2,p.SurroundSize,p.SurroundSize);
center_target = NaN(p.numTrials*2,p.CenterSize,p.CenterSize);
surround_target = NaN(p.numTrials*2,p.SurroundSize,p.SurroundSize);

for n = 1:p.numTrials
    center = (sin(p.Freq*2*pi/p.CenterSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.phase(n,1)));
    center_grating(n,:,:) = (center .* centerGaussian);
    
    surround = (sin(p.Freq_surround*2*pi/p.SurroundSize*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.phase(n,2)));
    surround_grating(n,:,:) = (surround .* surroundGaussian);
    
    targetcenter = (sin(p.Freq*2*pi/p.CenterSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.phase(n,3)));
    center_target(n,:,:) = (targetcenter .* centerGaussian);
    
    targetsurround = (sin(p.Freq_surround*2*pi/p.SurroundSize*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.phase(n,4)));
    surround_target(n,:,:) = (targetsurround .* surroundGaussian);
end

%%%WINDOW SETUP
[window,rect] = Screen('OpenWindow', Screens(2), p.Grey);
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
% load('MyGammaTable.mat');
% Screen('LoadNormalizedGammaTable', window, repmat(gammaTable, [1 3]));
HideCursor;
white = 255; green = [0 255 0];

% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


% Define coordinates where to draw the stimuli
CenterX = rect(3)/2; CenterY = rect(4)/2;
% coordinates for location on left and right side of fixation
Patch =  [CenterX CenterY];


Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);


%%% START THE EXPERIMENT
% Draw some text to the screen first outside of the experimental loop:

% experiment setup
space = zeros(1,256); space(KbName('Space')) = 1;
PsychHID('KbQueueCreate',deviceNumber, space);
PsychHID('KbQueueStart', deviceNumber);


WelcomeText = ['Welcome!' '\n' '\n'...
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
DrawFormattedText(window, WelcomeText, 'center', 'center', 255);
Screen('Flip', window);

while 1
    [pmbutton, ~] = PsychPowerMate('Get', powermate);
    if pmbutton == 1;
        break;
    end
end


Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
Screen('Flip', window);

StartTime = GetSecs;
WaitSecs(t.Starttime);

% make sure we can press esc to quit the experiment
StartKey=zeros(1,256); StartKey(KbName({'ESCAPE'})) = 1;
PsychHID('KbQueueCreate', deviceNumber, StartKey);


% Preallocate some variables
data.EstimatedContrast = NaN(p.numTrials, 1);
data.DifferenceContrast = NaN(p.numTrials, 1);
data.ResponseTime = NaN(p.numTrials, 1);
%%%%%TRIAL LOOP%%%%%%

for n = 1:p.numTrials
    
    CenterStimulus = Screen('MakeTexture', window, squeeze(center_grating(n,:,:))* (p.TrialEvents(n,2) * p.Grey) + p.Grey);
    
    surroundtext(:,:,1) = squeeze(surround_grating(n,:,:))* (p.TrialEvents(n,3)*p.Grey) + p.Grey;
    surroundtext(:,:,2) = transparencyMask;
    SurroundStimulus = Screen('MakeTexture', window, surroundtext);
    
    CheckerMaskpos = Screen('MakeTexture', window, FullChecker);
    CheckerMaskneg = Screen('MakeTexture', window, FullCheckerNeg);

    if p.TrialEvents(n,1) ~= 4 && p.TrialEvents(n,4) == 1
        Screen('DrawTexture', window, CenterStimulus, [], CenterRectOnPoint([0 0 p.CenterSize p.CenterSize], Patch(1), Patch(2)), p.TrialEvents(n,5))
    end
    
    if (p.TrialEvents(n,1) ~= 3 && p.TrialEvents(n,4) == 2) || (p.TrialEvents(n,1) == 2 && p.TrialEvents(n,4) == 2)
        Screen('DrawTexture', window, SurroundStimulus, [], CenterRectOnPoint([0 0 p.SurroundSize p.SurroundSize], Patch(1), Patch(2)), p.TrialEvents(n,5))
        Screen('FrameOval', window, p.Grey, CenterRectOnPoint([0 0 p.CenterSize+p.GapSize p.CenterSize+p.GapSize], Patch(1), Patch(2))', p.GapSize, p.GapSize)
    end
    
    Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation])
    Screen('Flip', window);
    WaitSecs(t.stimon);
    
   
       
    for f = 1: t.flickertime/(t.flicker*2)
            
        Screen('DrawTexture', window, CheckerMaskpos, [],CenterRectOnPoint([0 0 p.SurroundSize p.SurroundSize], Patch(1), Patch(2)));
        Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);

        Screen('DrawTexture', window, CheckerMaskneg, [],CenterRectOnPoint([0 0 p.SurroundSize p.SurroundSize], Patch(1), Patch(2)));
        Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);

    end
    
    % retention interval
    Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation])
    Screen('Flip', window);
    WaitSecs(t.retention);
    
    % Time to show stimulus in vwm condition    
    if p.TrialEvents(n,1) ~= 4 && p.TrialEvents(n,4) == 2
        Screen('DrawTexture', window, CenterStimulus, [], CenterRectOnPoint([0 0 p.CenterSize p.CenterSize], Patch(1), Patch(2)), p.TrialEvents(n,5))
    end
    
    if (p.TrialEvents(n,1) ~= 3 && p.TrialEvents(n,4) == 1) || (p.TrialEvents(n,1) == 2 && p.TrialEvents(n,4) == 1)
        Screen('DrawTexture', window, SurroundStimulus, [], CenterRectOnPoint([0 0 p.SurroundSize p.SurroundSize], Patch(1), Patch(2)), p.TrialEvents(n,5))
        Screen('FrameOval', window, p.Grey, CenterRectOnPoint([0 0 p.CenterSize+p.GapSize p.CenterSize+p.GapSize], Patch(1), Patch(2))', p.GapSize, p.GapSize)
    end
    Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation])    
    Screen('Flip', window);
    WaitSecs(t.stimon);
    
    for f = 1: t.flickertime/(t.flicker*2)
            
        Screen('DrawTexture', window, CheckerMaskpos, [],CenterRectOnPoint([0 0 p.SurroundSize p.SurroundSize], Patch(1), Patch(2)));
        Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);

        Screen('DrawTexture', window, CheckerMaskneg, [],CenterRectOnPoint([0 0 p.SurroundSize p.SurroundSize], Patch(1), Patch(2)));
        Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);

    end
    
    % retention interval
    Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation])
    Screen('Flip', window);
    WaitSecs(t.retention);    
    
    % set up button press
    PsychHID('KbQueueStart', deviceNumber);
    PsychHID('KbQueueFlush');
       
    % Show center or surround Grating and allow user to change contrast
    if p.TrialEvents(n,1) == 1 || p.TrialEvents(n,1) == 3
        Target = Screen('MakeTexture', window, squeeze(center_target(n,:,:))* (p.probecontrast(n)*p.Grey) + p.Grey);
        Screen('DrawTexture', window, Target, [], CenterRectOnPoint([0 0 p.CenterSize p.CenterSize], Patch(1), Patch(2)), p.TrialEvents(n,5))
    elseif p.TrialEvents(n,1) == 2 || p.TrialEvents(n,1) == 4
        Target = Screen('MakeTexture', window, squeeze(surround_target(n,:,:))* (p.probecontrast(n)*p.Grey) + p.Grey);
        Screen('DrawTexture', window, Target, [], CenterRectOnPoint([0 0 p.SurroundSize p.SurroundSize], Patch(1), Patch(2)), p.TrialEvents(n,5))
        Screen('FillOval', window, p.Grey, CenterRectOnPoint([0 0 p.CenterSize+p.GapSize p.CenterSize+p.GapSize], Patch(1), Patch(2))')
    end
    
    Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation])
    Screen('Flip', window);
    starttrial = GetSecs; % get the start time of each trial
    
    [~, startangle] = PsychPowerMate('Get', powermate);
    est_contrast = p.probecontrast(n);
    
    while 1 %start inf loop
        % Query PowerMate button state and rotation angle in "clicks"
        [pmbutton, angle] = PsychPowerMate('Get', powermate);
        % 1st button is the "or" of the 1st mouse button and the actual PowerMate button
        if startangle ~= angle
            
            % Convert turn of dial first to degrees and then to contrast:
            %             angles = (angle * 3.8298)/360;
            angles = ((startangle-angle)*3.8298);
            changecontrast = angles/360;
            est_contrast = est_contrast - changecontrast; % update the contrast relative to last dial position
            % Make sure we stay in range
            
            if est_contrast > 1
                est_contrast = 1;
            elseif est_contrast < 0
                est_contrast = 0.001;
            end
            if p.TrialEvents(n,1) == 1 || p.TrialEvents(n,1) == 3
                Target = Screen('MakeTexture', window, squeeze(center_target(n,:,:))* (est_contrast*p.Grey) + p.Grey);
                Screen('DrawTexture', window, Target, [], CenterRectOnPoint([0 0 p.CenterSize p.CenterSize], Patch(1), Patch(2)), p.TrialEvents(n,5))
            elseif p.TrialEvents(n,1) == 2 || p.TrialEvents(n,1) == 4
                Target = Screen('MakeTexture', window, squeeze(surround_target(n,:,:))* (est_contrast*p.Grey) + p.Grey);
                Screen('DrawTexture', window, Target, [], CenterRectOnPoint([0 0 p.SurroundSize p.SurroundSize], Patch(1), Patch(2)), p.TrialEvents(n,5))
                Screen('FillOval', window, p.Grey, CenterRectOnPoint([0 0 p.CenterSize+p.GapSize p.CenterSize+p.GapSize], Patch(1), Patch(2))')
            end
            Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation])
            Screen('Flip', window);
            
            startangle = angle;
        end
        if pmbutton == 1;
            data.EstimatedContrast(n) = est_contrast;
            if p.TrialEvents(n,1) == 1 || p.TrialEvents(n,1) == 3
                data.DifferenceContrast(n) = p.TrialEvents(n,2) - data.EstimatedContrast(n);
            else
                data.DifferenceContrast(n) = p.TrialEvents(n,3) - data.EstimatedContrast(n);
            end
            data.ResponseTime(n) = (GetSecs - starttrial);
            pmbutton = 0;
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
    
    %Present center fixation
    Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
    Screen('Flip', window);
    WaitSecs(t.iti);
    
    
    %%%Rest period
    if n== ceil(p.numTrials/2)
        rest = GetSecs;
        
        RestText = ['Halfway there! You can take a short break now, ' '' '\n' ...
            'or press the dial to continue' '\n' '\n' ];
        DrawFormattedText(window, RestText, 'center', 'center', white);
        Screen('Flip', window);
        pmbuttonbreak = 0;
        
        while 1
            [pmbuttonbreak, a] = PsychPowerMate('Get', powermate);
            if pmbuttonbreak == 1;
                break;
            end
        end
        Screen('FillOval', window, green, [CenterX-p.OuterFixation CenterY-p.OuterFixation CenterX+p.OuterFixation CenterY+p.OuterFixation]);
        Screen('Flip', window);
        WaitSecs(t.iti);
        t.Resttime = (GetSecs-rest)/60;
        
    end
    
end

t.EndTime = (GetSecs-StartTime)/60; %Get endtime of the experiment in seconds
%Draw some more text to the screen outside of the loop:
Screen(window,'TextSize',30);
ByebyeText = 'Great work! You have finished this run';
DrawFormattedText(window, ByebyeText, 'center', 'center', [255 255 255]);
Screen('Flip', window);
WaitSecs(2);
Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
Screen('CloseAll')

%%%SAVE OUT THE DATA FILE
cd(datadir);
TheData(runnumber).t = t;
TheData(runnumber).p = p;
TheData(runnumber).data = data;
eval(['save vWM_Comp_surrSuppression_', p.Subject, '.mat TheData'])

cd(expdir);
