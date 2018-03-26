function rms = calcRMS(x,winStep,winLength)
% x is the original signal
% winStep is the overlap between windows of length winLength
% RMS = sqrt(sum(F_i^2)/N) .. ie sqrt of mean of square of signal amplitude
% returns vector of RMS (power) values

nWindows=getNumWindows(x,winStep,winLength);

rms=zeros(nWindows,1);
%rms=zeros(length(x),1);

for i=1:nWindows
  xoffset = ((i-1)*winStep)+1;
  XX=x(xoffset:xoffset+winLength-1);
  % average and take the square root of each window of squared amplitude
  rms(i) = sqrt(mean(XX.^2));
%  rms(xoffset:i*winStep)=sqrt(mean(XX.^2)); % pad spaces with prev value
end

%% pad to end with last value
%if nWindows*winStep<length(x)
%  rms(nWindows*winStep+1:end)=rms(nWindows*winStep);
%end

