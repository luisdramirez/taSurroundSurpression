%%% createSurroundSurpressionStim
commandwindow
testingOn = 1; % option to test stimuli

if testingOn
    % Screen setup
    PsychDefaultSetup(2); % default values for PTB
    screenNumber = max(Screen('Screens')); % look at available screens
    
    % Defining white and grey
    white = WhiteIndex(screenNumber);
    grey = white / 2;
    
    % Skipping sync sests
    Screen('Preference', 'SkipSyncTests', 2);
    
    % Open the screen
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2, ...
        [], [], kPsychNeed32BPCFloat);
end

% --- 
% Grating Information
% ---

% Dimension of the region where the gabor will be in drawn (pixels)
gaborDimPix = windowRect(4) / 2;

% Sigma of Gaussian
sigma = gaborDimPix / 7;

% Obvious parameters
orientation = 0;
contrast = 0.8;
aspectRatio = 1.0;
phase = 0;

% Spatial frequency (cycles per pixel)
% One cycle = grey-black-grey-white-grey (one black and one white lobe)

numCycles = 5;
freq = numCycles / gaborDimPix;

backgroundOffset = [0.5 0.5 0.5 0];
disableNorm = 1;
preContrastMultiplier = 0.5;
gabortex = CreateProceduralGabor(window, gaborDimPix, gaborDimPix, [], ...
    backgroundOffset, disableNorm, preContrastMultiplier);

% Randomize phase of the Gabors and make a properties matrix
propertiesMat = [phase, freq, sigma, contrast, aspectRatio, 0, 0, 0];

% ---
% Draw gabor
% --- 

Screen('DrawTextures', window, gabortex, [], [], orientation, [], [], [], [], ...
    kPsychDontDoRotation, propertiesMat');

% Flip to screen
Screen('Flip', window);

% Wait for a button press to exit
KbWait;

% Clear screen
sca;


