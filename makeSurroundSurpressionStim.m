%%% makeTaSurroundSurpression
% 
function makeSurroundSurpressionStim(run)
%% Run setup
%run = ;
saveStim = 1;
saveFigs = 0;

%% Add paths

%% file i/o
stimDir = '/Users/luisramirez/Documents/Boston University/Ling Lab/taSurroundSurpression/stimuli'; % where to store stimulus file
stimFile = ['taSurroundSurpression' num2str(run)]; % name of stimulus file

%% Screen setup
Screens = Screen('Screens'); % look at available screens
p.ScreenWidthPixels = Screen('Rect', Screens(1));
ScreenWidth = 36; % 29 cm macbook air, 40 cm trinitron crt, 60 cm Qnix screen
ViewDistance = 68; % in cm, ideal distance: 1 cm equals 1 visual degree
VisAngle = (2*atan2(ScreenWidth/2, ViewDistance))*(180/pi); % Visual angle of the whole screen
p.ppd = round(p.ScreenWidthPixels(3)/VisAngle); % pixels per degree visual angle
p.Grey = 128;
%% Keys setup


%% Timing parameters

refrate = 60; % Hz
nFramesPerTarget = 8 ; 
targetDur = refrate/nFramesPerTarget; % (s)
targetSOA = 15/60; %18/60 = .300; 15/60 = .250; 16/60 = .267; %0.6; % (s)



%% Stimuli/Targets setup

%% Blocks setup



%% Sound setup
Fs = 44100;
cueFreqs = [1046.5 440]; % [higher high C = target 1, lower A = target 2]
for iTone = 1:numel(cueFreqs)
    tone = MakeBeep(cueFreqs(iTone), cueDur, Fs);
    cueTones(iTone,:) = applyEnvelope(tone, Fs);
end
%% Store all stimulus parameters

%% Make the stimuli

%% Setup targets

%% Create stimulus structure

end