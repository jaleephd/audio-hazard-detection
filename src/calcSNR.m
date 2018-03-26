function snr = calcsnr(x, n)
%
% calculate the signal to noise ratio using SNR = mean / SD
%
% params:
%   - x contains a vector of signal samples
%   - n is length of moving average window
% returns:
%   - snr the signal to noise ratio (a vector the same length as x)

[m,s]=calcMovingAvgStd(x,n);
snr=m./s;


