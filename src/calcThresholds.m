function thresholds = calcThresholds(yd,mWinLen,sdthresh)
%function thresholds = calcThresholds(yd,mWinLen,sdThresh,tsw) % for debugging
%
% yd is frequency components of signal over time, converted to dB:
%    - each row represents a time slice
%    - each column represents a frequency component
%    yd (time,freq) = amplitude of that frequency component in dB
% mWinLen, is the length of the moving window
% sdThresh is the threshold defined in terms of SD from mean
%
% returns a vector of threshold over time

% create a moving average and standard deviation of total signal over time,
% to look at spectral informaion that stands out, rather than just
% background noise, and use this to create the global threshold vector from
% smoothed mean dB of spectrum

sa=mean(yd'); % average dB across all freq
sd=std(yd'); % SD
b=ones(1,mWinLen) * 1/mWinLen;
ma=filter(b,1,sa'); % calculate moving average of mean
md=std(ma);

% create the threshold vector from smoothed mean dB of spectrum
thresholds=ma+(md*sdthresh);


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


