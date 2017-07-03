function [repsonse, timing] = runSurroundSuppression(run, subj)
%%

% INPUTS:
    % run number
    % subject name

% OUTPUTS:
    % response
    % timing 

%% Window Setup
[window,rect] = Screen('OpenWindow', Screens(1), p.Grey);
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
load('MyGammaTable.mat');
Screen('LoadNormalizedGammaTable', window, repmat(gammaTable, [1 3]));
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

%% Run Experiment


end