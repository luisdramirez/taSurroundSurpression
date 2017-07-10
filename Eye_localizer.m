
%%% Localizer 

% Script shows flickering checkerboard, 
% Subjects perform a RSVP task
% Last adjusted 1/27/15 - Ilona Bloem

% 136 TRs

% RSVP task - index finger for 'j', middle finger for 'k'

% try
%%% PREPARE AND COLLECT INFO
echo off
clear all
KbName('UnifyKeyNames')
Screen('Preference', 'SkipSyncTests', 1);
input('hit enter to begin...  ');

%%%%%%%%

output.Subject = 'tmp';

%%%%%%%%

realScan = 0; % 1 = yes, 0 = no

%%%%%%%%

if realScan == 0
    eyeTrackON = 0; % 1 == start eye tracking, 0 = no eyetracking
else
    eyeTrackON = 1; % 1 == start eye tracking, 0 = no eyetracking
end

% scriptpath = '/Users/thelinglab/Desktop/Ilona/Tuned_Suppression';
% cd(scriptpath);

t.MySeed = sum(100*clock);
rand('state', t.MySeed); %rng(t.MySeed);

[keyboardIndices, productNames, ~] = GetKeyboardIndices;
if realScan == 1
    deviceString= 'Teensy Keyboard/Mouse';
else
    deviceString = 'Apple Internal Keyboard / Trackpad'
end
for i=1:length(productNames)%for each possible device   
    if strcmp(productNames{i},deviceString)%compare the name to the name you want      
        deviceNumber=keyboardIndices(i);%grab the correct id, and exit loop       
        break;        
    end    
end
if deviceNumber==0%%error checking
    error('No device by that name was detected');    
end
triggerKey = 46;                    % KbName('=+'); 
keylist = ones(1,256);              % keys for KbQueueCreate
keylist(triggerKey) = 0;            % dont want to record trigger
if realScan == 0
    keyPressNumbers = {'80', '79'}; % {'80', '79'} for arrows on macbook
else
    keyPressNumbers = {'30', '31'}; % for scanner, 
end
   

%%% SCREEN PARAMETERS
w.whichScreen = 0;
if realScan == 0
    w.ScreenWidth = 29;             % horizontal display size - 41.5 in scanner;
    w.ViewDistance = 57;            % in cm, ideal distance: 1 cm equals 1 visual degree (at 57 cm) - 107.5 at scanner with eye-tracking, 98 normal screen
else
    w.ScreenWidth = 41.5;           % horizontal display size - 41.5 in scanner, 44.5 without eye-tracking;
    w.ViewDistance = 107.5;         % in cm, ideal distance: 1 cm equals 1 visual degree (at 57 cm) - 107.5 at scanner with eye-tracking, 88 normal screen
end
w.frameRate = 60;
w.ScreenSizePixels = Screen('Rect', w.whichScreen); %Scanner display = [0 0 1024 768]; 
w.VisAngle = (2*atan2(w.ScreenWidth/2, w.ViewDistance))*(180/pi); % Visual angle of the whole screen
stim.ppd = round(w.ScreenSizePixels(3)/w.VisAngle); % pixels per degree visual angle
sizepixel = (2*atan2((w.ScreenWidth/w.ScreenSizePixels(3))/2, w.ViewDistance))*(180/pi); % in visual degree
fNyquist = 0.5/sizepixel;

%%% TIMING PARAMETERS
t.TheDate = datestr(now,'yymmdd');  %Collect todays date
t.TimeStamp = datestr(now,'HHMM');  %Timestamp for saving out a uniquely named datafile (so you will never accidentally overwrite stuff)
t.TR = 2;                           % TR length
t.RSVP_duration = .2;               % letter presentation duration
t.initon = 0;                       % 14 s baseline
t.finalon = 0;
t.IBI = 0;                          % Inter Block Interval
t.stimon = .25;                     % grating presentation duration
t.stimoff = .25;
t.on_dur = 16;                      % on block duration  (s)
t.responseDur = 1.5;                % response time
t.Repetitions = 6;                  % number of repetitions of each condition
t.numblocks = t.Repetitions*2;      % number of blocks (2 attention)
t.nPresentations = (t.on_dur-t.IBI)/(t.stimon+t.stimoff); 

output.TotalTRs = (t.initon + t.on_dur * (t.numblocks+1)) / t.TR;
output.TotalSecs = output.TotalTRs*t.TR;
stim.presentationsPerBlock = (t.on_dur-t.IBI)/t.stimon;
RefreshDur = 1/w.frameRate;
NumStimOnFrames = t.stimon/RefreshDur;
NumStimOffFrames = t.stimoff/RefreshDur;
WhenStart = t.initon + t.IBI : t.on_dur : (t.numblocks+1)*(t.on_dur)+t.initon;
WhenStop = WhenStart+(t.on_dur-t.IBI);

%%% STIMULUS PARAMETERS
stim.size = 15;         % in visual degree
stim.annulus = 3;
stim.fixation = 0.7;
stim.diameter = round(stim.size*stim.ppd); % in pixels
stim.annulus_diam = round(stim.annulus*stim.ppd);
stim.fix_diam = round(stim.fixation*stim.ppd);
stim.outer_fixation = round(1.4 * stim.fix_diam);

% stim.orientations = [45 135]; % two obliques
stim.Freq = 0.5;
stim.Contrast = 1; % max is 50% -> adding two signal together
stim.Grey = 127;
stim.Ampl = stim.Contrast * stim.Grey;

stim.distractorLetters = ['X' 'L' 'V' 'H' 'S' 'A' 'C' 'P' 'Z' 'Y'];
stim.TargetLetter = ['J', 'K'];
stim.p_target = 0.3;
Freq = stim.Freq .* stim.size;

%%% INITIATE & GENERATE DATA FILE
if exist(['Localizer_', output.Subject, '.mat']);
    load (['Localizer_', output.Subject, '.mat']);
    runnumber = length(TheData)+1;
    output.observer = [output.Subject '_L0' num2str(runnumber)];
else
    runnumber = 1;
    output.observer = [output.Subject '_L0' num2str(runnumber)];
end

%%% CREATE CHECKERBOARDS
% Gaussian Mask
xysize = round(stim.diameter);
[x,y] = meshgrid((-xysize/2):(xysize/2)-1, (-xysize/2):(xysize/2)-1);
gaussian_std = round(stim.ppd*2);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
Gaussian = zeros(xysize); Gaussian(eccen <= (xysize/2)-stim.ppd/2) = 1; 
% Gaussian = conv2(Gaussian, fspecial('gaussian', stim.ppd, stim.ppd), 'same');

% Gaussian Annulus
[X,Y] = meshgrid(1:xysize,1:xysize);
Annulus = zeros(xysize);
r_eccen = sqrt((X-xysize/2).^2+(Y-xysize/2).^2); 	% calculate eccentricity of each point in grid relative to center
Annulus(r_eccen > stim.annulus_diam/2) = 1;
% Annulus = conv2(Annulus, fspecial('gaussian', 30, 10), 'same');
Annulus(r_eccen > stim.annulus_diam) = 1;
 
% Checkerboards
Mask = Gaussian .* Annulus;


radius_concentric = (sqrt(x.^2+y.^2));
concentric = log(radius_concentric) * 30*pi/180;
concentric_radial = log(radius_concentric) * 0*pi/180;
radial_1 = 1 * (atan2(x,y)) + 0;
radial_2 = -1 * (atan2(x,y)) + 0;
grating_1 = (cos((concentric+radial_1) * Freq));
grating_2 = (cos((concentric+radial_2) * Freq));
grating_3 = (cos((concentric_radial+radial_2) * Freq*12));
grating = grating_1+grating_2+grating_3;
checks_tmp = sign(grating./max(max(grating)));

% gratingUp = (square(Freq*2*pi/stim.diameter*(X.*sin(pi/2)+Y.*cos(pi/2))));
% gratingDown = (square(Freq*2*pi/stim.diameter*(X.*sin(pi)+Y.*cos(pi))));
% CheckerBoard1 = (gratingUp .* gratingDown) .* Mask;
CheckerBoard1 = checks_tmp .* Mask;
CheckerBoard2 = (-1*CheckerBoard1);
CheckerBoard1 = CheckerBoard1*stim.Ampl+stim.Grey;
CheckerBoard2 = CheckerBoard2*stim.Ampl+stim.Grey;

%Trialevents
stim.TrialEvents = [repmat([0; 1], t.Repetitions, 1); 0];

%%% WINDOW SETUP
AssertOpenGL;
[window rect] = Screen('OpenWindow',w.whichScreen, stim.Grey);
HideCursor;

%%% MAKE COLOR LOOKUP TABLE AND APPLY GAMMA CORRECTION
OriginalCLUT = Screen('LoadCLUT', window);
yellow = [255 255 0]; red = [255 0 0]; green = [0 255 0]; blue = [0 0 255];
white = WhiteIndex(w.whichScreen);
black = BlackIndex(w.whichScreen);
if realScan == 1
    MyCLUT = load('linearizedCLUT.mat');
    Screen('LoadNormalizedGammaTable', window, MyCLUT.linearizedCLUT);
end

%%%%%%%%%%%%%%%%%%%
%%% If eyeTrackOn == 1, eye link setup
if eyeTrackON == 1
    [el edf_filename] = eyeTrackingOn(window, output.observer, rect, stim.ppd);
end
%%%%%%%%%%%%%%%%%%%

%%% CREATE PATCHES
CenterX = w.ScreenSizePixels(3)/2; CenterY = w.ScreenSizePixels(4)/2;
CenterPatch = [CenterX-stim.diameter CenterY-stim.diameter CenterX+stim.diameter CenterY+stim.diameter];
Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);
bbox = Screen('TextBounds', window, 'X');
newRect = CenterRectOnPoint(bbox, CenterX, CenterY);
tx = newRect(1); ty = newRect(2);

% Build designmatrix - 
% 1:baseline; 
% 2:stimon
output.designMatrix = zeros(2,output.TotalTRs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% START THE EXPERIMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialize run
tCount = 1;
block_counter = 0;
RSVP_counter = 0;
RSVP_time = 0;
WhatLetter = RandSample(stim.distractorLetters);
RSVP_rightwrong = 0;
responsewindowTrig_RSVP = 0;

% wait for backtick sync with scanner
Screen('FillOval', window, white, CenterRectOnPoint([0 0 stim.outer_fixation stim.outer_fixation], CenterX, CenterY));
Screen('FillOval', window, stim.Grey, CenterRectOnPoint([0 0 stim.fix_diam stim.fix_diam], CenterX, CenterY));
Screen('DrawText', window, '~', tx, ty + 40, 255, 0);
Screen('Flip', window);
% WaitSecs(2);
% GetClicks;

KbTriggerWait(triggerKey, deviceNumber);

PsychHID('KbQueueCreate', deviceNumber, keylist);

%%%%%%%%%%%%%%%%%%%
%%% If eyeTrackOn == 1, start recording
if eyeTrackON == 1
    [status el] = eyeTrackingRecord(el, rect, stim.ppd);
end
%%%%%%%%%%%%%%%%%%%

%%% INTIAL BASELINE----------------------------------------------------
StartTime = GetSecs;
tic;
Screen('FillOval', window, white, CenterRectOnPoint([0 0 stim.outer_fixation stim.outer_fixation], CenterX, CenterY));
Screen('FillOval', window, stim.Grey, CenterRectOnPoint([0 0 stim.fix_diam stim.fix_diam], CenterX, CenterY));
Screen('Flip', window, StartTime+WhenStart(1)-t.IBI, [], [], 1);

output.trueBaseline = toc;
temp(tCount) = GetSecs-StartTime;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STIMULI LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:t.numblocks+1
    allowtarget = 0;
    responsewindowTrig = 0;
    responsewindowTrig_RSVP = 0;
    block_counter = block_counter+1;
    Rotate_it = 0;
    if rand > 0.5, rotate = 5; else rotate = -5; end % randomise which way to rotate
    
    Screen('Flip', window, [], [], 1);

    if stim.TrialEvents(ii,1) == 1 % stimulus on
        baseline = 0;
        output.designMatrix(2,(ii-1)*t.on_dur/t.TR+1:(ii)*t.on_dur/t.TR) = ones(1,t.on_dur/t.TR);
    else % blank period
        baseline = 1; 
        output.designMatrix(1,(ii-1)*t.on_dur/t.TR+1:(ii)*t.on_dur/t.TR) = ones(1,t.on_dur/t.TR);
    end
    Screen('Flip', window, StartTime+WhenStart(block_counter), [], [], 1);
    
    % save real timing
    tCount = tCount+1;
    output.trueOn_Time(tCount) = GetSecs - StartTime;

    tic;
    for n = 1:t.nPresentations 
        Rotate_it = Rotate_it + rotate;       
        
        if baseline == 0 % not a baseline block
            Checker1 = Screen('MakeTexture', window, CheckerBoard1);
            Checker2 = Screen('MakeTexture', window, CheckerBoard2);
        end      
        
        % tic
        for kk = 1:NumStimOnFrames
            if GetSecs < StartTime+WhenStart(block_counter)+1 || GetSecs > StartTime+WhenStop(block_counter)-2
                allowtarget = 0;
            end
                       
            if GetSecs > StartTime+WhenStop(block_counter)
                break
            end
            
            if baseline == 0
                Screen('DrawTexture', window, Checker1, [], [], Rotate_it);                    
            end
            Screen('FillOval', window, white, CenterRectOnPoint([0 0 stim.outer_fixation stim.outer_fixation], CenterX, CenterY));
            Screen('FillOval', window, stim.Grey, CenterRectOnPoint([0 0 stim.fix_diam stim.fix_diam], CenterX, CenterY));

            if GetSecs > RSVP_time
                allowtarget = allowtarget+1;
                if allowtarget > 10 %&& baseline == 0
                    RSVPTargetOrNot = rand < stim.p_target;
                else
                    RSVPTargetOrNot = 0;
                end
                if RSVPTargetOrNot
                    letter_indx = RandSample([1 2]);
                    WhatLetter = stim.TargetLetter(letter_indx);
                    allowtarget = 0;
%                     if stim.TrialEvents(ii,1) == 1                  
                        PsychHID('KbQueueStart', deviceNumber);
                        responsewindowTrig_RSVP = 1;
                        responsewindow_RSVP = GetSecs+t.responseDur;
                        RSVP_counter = RSVP_counter+1;
                        RSVP_rightwrong(RSVP_counter) = 0;
                        output.whatletter(RSVP_counter) = WhatLetter;
%                     end
                    
                end
                if ~RSVPTargetOrNot
                    WhatLetter_tmp = RandSample(stim.distractorLetters);
                    while strcmp(WhatLetter, WhatLetter_tmp)
                        WhatLetter_tmp = RandSample(stim.distractorLetters);
                    end
                    WhatLetter = WhatLetter_tmp;
                end
                RSVP_time = GetSecs + t.RSVP_duration;
            end
            Screen('DrawText', window, WhatLetter, tx, ty, 255, 0);
            Screen('Flip', window, 0, [], 1);
                    
            if responsewindowTrig_RSVP == 1;
                if GetSecs >= responsewindow_RSVP
                    [pressed, firstpress] = PsychHID('KbQueueCheck', deviceNumber);
                    whichkeys = find(firstpress);
                    if ((strcmp(num2str(whichkeys), keyPressNumbers{1}) && letter_indx==1) || (strcmp(num2str(whichkeys), keyPressNumbers{2}) && letter_indx==2))
                        RSVP_rightwrong(RSVP_counter) = 1
                    else
                        RSVP_rightwrong(RSVP_counter) = 0;
                    end
                    KbQueueFlush();
                    responsewindowTrig_RSVP = 0;
                end
            end
        end
        % onframes(n) = toc; 
        
%         Rotate_it = Rotate_it +5;       
        
        % tic
        for kk = 1:NumStimOffFrames   
            
            if GetSecs < StartTime+WhenStart(block_counter)+1 ||  GetSecs > StartTime+WhenStop(block_counter)-2
                allowtarget = 0;
                allowtarget_orient = 0;
            end
            
            if GetSecs > StartTime+WhenStop(block_counter)
                break
            end
            
            if GetSecs > RSVP_time
                allowtarget = allowtarget+1;
                if allowtarget > 10 %&& baseline == 0
                    RSVPTargetOrNot = rand < stim.p_target;
                else
                    RSVPTargetOrNot = 0;
                end
                if RSVPTargetOrNot
                    letter_indx = RandSample([1 2]);
                    WhatLetter = stim.TargetLetter(letter_indx);
                    allowtarget = 0;
%                     if stim.TrialEvents(ii,1) == 1
                        PsychHID('KbQueueStart', deviceNumber);
                        responsewindowTrig_RSVP = 1;
                        responsewindow_RSVP = GetSecs+t.responseDur;
                        RSVP_counter = RSVP_counter+1;
                        RSVP_rightwrong(RSVP_counter) = 0;
                        output.whatletter(RSVP_counter) = WhatLetter;
%                     end
                    
                end
                if ~RSVPTargetOrNot
                    WhatLetter_tmp = RandSample(stim.distractorLetters);
                    while strcmp(WhatLetter, WhatLetter_tmp)
                        WhatLetter_tmp = RandSample(stim.distractorLetters);
                    end
                    WhatLetter = WhatLetter_tmp;
                end
                RSVP_time = GetSecs + t.RSVP_duration;
            end
            
            if baseline == 0
                Screen('DrawTexture', window, Checker2, [], [], Rotate_it);                    
            end
            
            Screen('FillOval', window, white, CenterRectOnPoint([0 0 stim.outer_fixation stim.outer_fixation], CenterX, CenterY));
            Screen('FillOval', window, stim.Grey, CenterRectOnPoint([0 0 stim.fix_diam stim.fix_diam], CenterX, CenterY));            
            Screen('DrawText', window, WhatLetter, tx, ty, 255, 0);
            Screen('Flip', window, 0, [], 1);
            
            if responsewindowTrig_RSVP == 1;
                if GetSecs >= responsewindow_RSVP
                    [pressed, firstpress] = PsychHID('KbQueueCheck', deviceNumber);
                    whichkeys = find(firstpress);
                    if ((strcmp(num2str(whichkeys), keyPressNumbers{1}) && letter_indx==1) || (strcmp(num2str(whichkeys), keyPressNumbers{2}) && letter_indx==2))
                        RSVP_rightwrong(RSVP_counter) = 1
                    else
                        RSVP_rightwrong(RSVP_counter) = 0;
                    end
                    KbQueueFlush();
                    responsewindowTrig_RSVP = 0;
                end
            end
            
        end    
        % offframes(n) = toc; 
        output.trueOn_duration(block_counter) = toc;
        
    end
    
    if baseline == 0;
        Screen('Close',Checker1);
        Screen('Close',Checker2);
    end

end

output.totalScanDuration = GetSecs - StartTime;

%%% make & output
output.RSVP_rightwrong = RSVP_rightwrong;
output.mean_RSVP_rightwrong = mean(RSVP_rightwrong);
TheData{runnumber}.stim = stim;
TheData{runnumber}.t = t;
TheData{runnumber}.w = w;
TheData{runnumber}.mylog.stimtimes_s = {output.trueOn_Time};
TheData{runnumber}.mylog.durationss_s = {output.trueOn_duration};
TheData{runnumber}.mylog.designMatrix = output.designMatrix;
TheData{runnumber}.output = output;
eval(['save ', 'Localizer_', output.Subject, '.mat TheData']);

Screen('FillOval', window, white, CenterRectOnPoint([0 0 stim.outer_fixation stim.outer_fixation], CenterX, CenterY));
Screen('FillOval', window, stim.Grey, CenterRectOnPoint([0 0 stim.fix_diam stim.fix_diam], CenterX, CenterY));            
Screen('DrawText', window, 'DONE.', CenterX-20, CenterY-stim.annulus_diam/2);
Screen('Flip', window, 0, [], 1);

WaitSecs(3)

%%% close screen
ShowCursor
Screen('CloseAll')


%%%%%%%%%%%%%%%%%%%
if eyeTrackON == 1
    Eyelink('StopRecording');
    Eyelink('CloseFile');
    Eyelink('ReceiveFile',edf_filename);
end
%%%%%%%%%%%%%%%%%%%

disp('------------------------------------------------');
disp(['RSVP Accuracy: ' num2str(mean(RSVP_rightwrong)), ' (', num2str(sum(RSVP_rightwrong)), '/', num2str(numel(RSVP_rightwrong)), ')']);
disp('------------------------------------------------');


% end


