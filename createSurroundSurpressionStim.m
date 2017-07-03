%%% createSurroundSurpressionStim
commandwindow
testingOn = 1; % option to test stimuli

if testingOn
    % Screen setup
    PsychDefaultSetup(2); % default values for PTB
    screenNumber = max(Screen('Screens')); % fix display to most external display
    param.screenWidthPixels = Screen('Rect', screenNumber); % screen width in pixels 
    screenWidth = [];
    viewDistance = [];
    visAngle = (2*atan2(screenWidth/2, viewDistance))*(180/pi); % visual angle of the whole screen
    param.pixPerDeg = round(p.screenWidthPixels(3)/visAngle); % pixels per degree visual angle
    
    % Defining colors
    param.white = WhiteIndex(screenNumber);
    param.grey = white / 2;
    param.green = [0 255 0];
    
    
    % Skipping sync sests
    %Screen('Preference', 'SkipSyncTests', 2);
    
    % Open the screen
    [window, windowRect] = Screen('OpenWindow', screenNumber, param.grey);
    originalCLUT = Screen('ReadNormalizedGammeTable', window);
    % load('MyGammaTable.mat')
    % Screen('LoadNormalizedGammaTable', window, repmat(gammaTable, [1
    % 3]));
    
    HideCursor;
    
    % Enable alpha blending
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Define where to draw stimuli
    centerX = rect(3)/2; centerY = rect(4)/2;
    
    % Coordinates for location on left and right side of fixation
    patch = [centerX centerY];
    
    Screen('TextStyle', window, 1);
    Screen('TextSize', window, 16);
end


%% timing
timing.stimOn = 1; % stimulus on screen (s)
timing.iti = 0.3; % inter-trial interval(s)
timing.flickerTime = 0.2; % (s)
timing.flicker = 0.025; % (s)
%% grating parameters
param.stimConfigurations = 1; % 1 = simultaneous, 2 = sequentially - center first, 3 = sequentially - surround first
param.numContrasts = 1;
param.testContrasts = 10.^linspace(log10(0.1),log10(0.75), param.numContrasts);
param.surroundContrasts = param.testContrasts; % surround and center should have same contrast

% size parameters
param.centerSize = round(1 * param.pixPerDeg);
param.surroundSize = param.centerSize * 3;
param.gapSize = round(0.08 * param.pixPerDeg); % space between inner and outer stimuli
param.outerFixation = round(0.05 * param.pixPerDeg);
param.innerFixation = param.outerFixation / 1.5;

freq = 2;
param.freq = param.centerSize / param.pixPerDeg * freq;
param.freqSurround = param.surroundSize / param.pixPerDeg * freq;
param.orientation = 0;
param.phase = randsample(1:180, param.numTrials*4, true);
param.phase = reshape(param.phase, [param.numTrials 4]);
param.probecontrast = randsample(0.1:0.01:0.9, param.numTrials, true);
param.orientationChecker = [0 90];
param.phaseChecker = [0 180];

%% create stimuli

% make mask to create circle for the center grating
[x,y] = meshgrid((-param.centerSize/2):(param.centerSize/2)-1, (-param.centerSize/2):(param.centerSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
centerGaussian = zeros(param.centerSize); centerGaussian(eccen <= (param.centerSize/2)) = 1;
% Gaussian = conv2(Gaussian, fspecial('gaussian', p.ppd, p.ppd), 'same');

% make mask to create circle for the surround grating
[x,y] = meshgrid((-param.surroundSize/2):(param.surroundSize/2)-1, (-param.surroundSize/2):(param.surroundSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
surroundGaussian = zeros(param.surroundSize); surroundGaussian(eccen <= (param.surroundSize/2)) = 1;

% Make transparency mask for aplha blending the two images
transparencyMask = zeros(param.surroundSize); transparencyMask(eccen >= ((param.centerSize)/2)) = 255;

% make unique grating for every trial
[Xc,Yc] = meshgrid(0:(param.centerSize-1),0:(param.centerSize-1));
[Xs,Ys] = meshgrid(0:(param.surroundSize-1),0:(param.surroundSize-1));

% Make checkerboard
checker1 = (square(p.freqSurround*2*pi/param.surroundSize*(Xs.*sin(param.orientationChecker(1)*(pi/180))+Ys.*cos(param.orientationChecker(1)*(pi/180)))-param.phaseChecker(1)));
checker2 = (square(p.freqSurround*2*pi/param.surroundSize*(Xs.*sin(param.orientationChecker(2)*(pi/180))+Ys.*cos(param.orientationChecker(2)*(pi/180)))-param.phaseChecker(2)));
fullChecker = (checker1.*checker2) .* surroundGaussian;
fullCheckerNeg = (fullChecker*-1);
fullChecker = fullChecker * (param.grey-1) + param.grey;
fullCheckerNeg = fullCheckerNeg * (param.grey-1) + param.grey;

% Make actual gratings
centerGrating = NaN(param.numTrials*2, param.centerSize, param.centerSize);
surroundGrating = NaN(param.numTrials*2, param.surroundSize,param.surroundSize);
centerTarget = NaN(param.numTrials*2, param.centerSize, param.centerSize);
surroundTarget = NaN(param.numTrials*2, param.surroundSize, param.surroundSize);

for n = 1:param.numTrials
    center = (sin(param.freq*2*pi/param.centerSize*(Xc.*sin(param.orientation*(pi/180))+Yc.*cos(param.orientation*(pi/180)))-param.phase(n,1)));
    centerGrating(n,:,:) = (center .* centerGaussian);
    
    surround = (sin(param.freqSurround*2*pi/param.surroundSize*(Xs.*sin(param.orientation*(pi/180))+Ys.*cos(param.orientation*(pi/180)))-param.phase(n,2)));
    surroundGrating(n,:,:) = (surround .* surroundGaussian);
    
    targetCenter = (sin(param.freq*2*pi/param.centerSize*(Xc.*sin(param.orientation*(pi/180))+Yc.*cos(param.orientation*(pi/180)))-param.phase(n,3)));
    centerTarget(n,:,:) = (targetcenter .* centerGaussian);
    
    targetSurround = (sin(param.freqSurround*2*pi/param.surroundSize*(Xs.*sin(param.orientation*(pi/180))+Ys.*cos(param.orientation*(pi/180)))-param.phase(n,4)));
    surroundTarget(n,:,:) = (targetSurround .* surroundGaussian);
end



