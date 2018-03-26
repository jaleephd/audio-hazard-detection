function scaled = scaleSignalToSD(x,mavWinLen,msdWinLen)
%
% scale signal in terms of standard deviation from moving average over time,
% to look at information that stands out, rather than background noise
%
% params:
% - x contains a vector of signal samples
% - mavWinLen, is the length of the moving window for average
% - msdWinLen, is the length of the moving window for std deviation
%
% returns
% - scaled, a vector of the signal scaled in terms of SD from moving average


[m,junk]=calcMovingAvgStd(x,mavWinLen);

% this determines whether to use moving SD (emulate real time processing)
% if data is only over a limited time frame (ie within a scenario)
% then can use SD over all time
if msdWinLen>0
  [junk,s]=calcMovingAvgStd(x,msdWinLen);
else
  s=std(x);
end


x2=x-m; % shift values relative to moving average to look for values +/- SD
scaled=x2./s; % scale in terms of SD

