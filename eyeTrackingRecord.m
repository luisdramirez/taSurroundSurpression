
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      EYE LINK SETUP                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [status el] = eyeTrackingRecord(el, rect, pixelPerDeg)

center = [rect(3)/2 rect(4)/2];
box = 2.5*pixelPerDeg;


 %%%%eyetracker stuff==============================================
% Send trial id message to Eyelink file

% Sending a 'TRIALID' message to mark the start of a trial in Data 
% Viewer.  This is different than the start of recording message 
% START that is logged when the trial recording begins. The viewer
% will not parse any messages, events, or samples, that exist in 
% the data file prior to this message. 
Eyelink('Message', 'TRIALID %d', 1);

% This supplies the title at the bottom of the eyetracker display
Eyelink('command', 'record_status_message "TRIAL %d"', 1); 
% Before recording, we place reference graphics on the host display
% Must be offline to draw to EyeLink screen
Eyelink('Command', 'set_idle_mode');
% clear tracker display and draw box at fix point
Eyelink('Command', 'clear_screen 0')
% Eyelink('command', 'draw_box %d %d %d %d 15', width/2-50, height/2-50, width/2+50, height/2+50);
Eyelink('command', 'draw_box %d %d %d %d 15',round(center(1)-box), round(center(2)-box), round(center(1)+box), round(center(2)+box));

Eyelink('Command', 'set_idle_mode');
WaitSecs(0.01);%.05


[status]=Eyelink('StartRecording');
Eyelink('Message', 'SYNCTIME');