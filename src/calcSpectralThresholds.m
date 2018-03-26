function thresholds = calcSpectralThresholds(yd,mavWinLen,msdWinLen,sdthresh)
%function thresholds = calcSpectralThresholds(yd,mavWinLen,msdWinLen,sdThresh,tsw) % for debugging
%
% yd is frequency components of signal over time, converted to dB:
%    - each row represents a time slice
%    - each column represents a frequency component
%    yd (time,freq) = amplitude of that frequency component in dB
% mavWinLen, is the length of the moving window for average
% msdWinLen, is the length of the moving window for std dev
% sdThresh is the threshold defined in terms of SD from mean
%
% returns a vector of threshold over time

% create a moving average and standard deviation of total signal over time,
% to look at spectral informaion that stands out, rather than just
% background noise, and use this to create the global threshold vector from
% smoothed mean dB of spectrum

sa=mean(yd'); % average dB across all freq
sd=std(yd'); % SD across freq components

%b=ones(1,mWinLen) * 1/mWinLen;
%ma=filter(b,1,sa'); % calculate moving average of mean

% calculate moving average and SD of mean
% Note: SD can be over all time or varying over time (based on windowed samples)
%       for the latter, use the value returned from calcMovingAvgStd()
[ma,junk] = calcMovingAvgStd(sa',mavWinLen);
if msdWinLen>0
  [junk,md] = calcMovingAvgStd(sa',msdWinLen);
  % create the threshold vector from smoothed mean dB of spectrum
  thresholds=ma+(md.*sdthresh);
else
  % if data is over a limited time frame (ie within a scenario)
  % then use SD over all time
  md=std(ma); % SD of moving average over all time
  % create the threshold vector from smoothed mean dB of spectrum
  thresholds=ma+(md*sdthresh);
end



if 0
clf reset % clear graph
hold on
title('Mean & SD Signal Amplitude');
xlabel('time (sec)');
ylabel('Amplitude');
plot(tsw(:),sa,'k-');
plot(tsw(:),ma,'b-');
lgn=legend('instantaneous','smoothed',4);

%plot(tsw(:),sa-sd,'g:');
%plot(tsw(:),sa+sd,'g:');
plot(tsw(:),sa-(sd*sdthresh),'r');
plot(tsw(:),sa+(sd*sdthresh),'r');
%plot(tsw(:),ma-md,'c:');
%plot(tsw(:),ma+md,'c:');
plot(tsw(:),ma-(md*sdthresh),'m');
plot(tsw(:),ma+(md*sdthresh),'m');
grid on;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
end


