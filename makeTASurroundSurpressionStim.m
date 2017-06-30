%%% makeTaSurroundSurpression

function makeTASurroundSurpressionStim(run)
%% Run setup
%run = ;
saveStim = 1;
saveFigs = 0;

%% Add paths

%% file i/o
stimDir = '/Users/luisramirez/Documents/Boston University/Ling Lab/taSurroundSurpression/stimuli';
stimFile = ['taSurroundSurpression' num2str(run)]; %name of stimulus file

%% Screen setup

%% Keys setup


%% Timing parameters

refrate = 60; % Hz
nFramesPerTarget = 8 ; 
targetDur = refrate/nFramesPerTarget; % (s)
targetSOA = 15/60; %18/60 = .300; 15/60 = .250; 16/60 = .267; %0.6; % (s)



%% Target setup

%% Blocks setup

%% Stimulus setup

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