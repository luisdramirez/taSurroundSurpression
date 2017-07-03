%%% loadDisplayParams

function [display] = loadDisplayParams(displayName)

switch displayName
    case 'luisMacbook'
        display.numPixels = 
        display.pixelSize = 
        display.distance = 
    otherwise
        error('Display name not recognized.')
end