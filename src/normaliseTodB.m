function [Yd,Yn,maxampl,minampl] = normaliseToDb(Y,maxampl,ceiling,floordB)
% Y is signal amplitudes from FFT, maxampl is previous maximum, celing is a hard amplitude ceiling (-1 if not using), floordB is the minimum allowed dB (eg -100)

if ceiling<0 % if not using a hard ceiling
  maximum = max(max(abs(Y))); % find the largest amplitude
  if maximum>maxampl     % set the maximum amplitude to be the highest
    maxampl=maximum;     % of the supplied previous maximum, or this max
  end
else
  maxampl=ceiling;
end

minampl=min(min(abs(Y))); % this should be (close to 0)
mratio=maxampl-minampl;

if maxampl~=0
  Yn=(abs(Y)-minampl)/mratio; % normalise to 0 - 1
else
  Yn=abs(Y);
end


% convert amplitude to dB
% taken from
%   http://hans.fugal.net/blog/2009/04/02/dft-magnitudepower-spectra
%  and
%   myspecgram.m by Paul Kienzle; modified by Sean Fulop March 2002

Yd = 20*log10(Yn);

% clip everything below floordB dB (get rid of -Inf). floor dB should be < 0
Yd=max(Yd,floordB);


