function res = showABD(ABDvar)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if ABDvar.titration_BCcurve
   res = plotABDarea(ABDvar,1);
else
   res = plotABDarea(ABDvar,0);
end
end