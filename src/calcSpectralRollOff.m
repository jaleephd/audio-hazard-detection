function [rolloff,rolloffidx] = calcSpectralRollOff(yn,freqs)

% yn (time,freq) is normalised amplitudes of the frequency component (0-1)
% freqs is a list of frequencies in the spectrum (or band)
%
% the spectral roll off frequency is where 95% of the signal energy
% (summed squared amplitude) is contained below this frequency
% ie
%   sum_0,fr a^2(f) = 0.95 sum_0_sr/2 a^2(f)
%   where fr is roll off frequency, sr/2 is the Nyquist frequency
%                                   (1/2 the sampling rate)

threshold=0.95;

rolloff=zeros(size(yn,1),1); % for each time slice
rolloffidx=zeros(size(yn,1),1); % for each time slice
for i=1:size(yn,1)
  total=sum(yn(i,:).^2);
  cuml=cumsum(yn(i,:).^2);
  ridx=min(find((cuml>=total*threshold))); % find threshold point
  rolloff(i)=freqs(ridx);
  rolloffidx(i)=ridx;
end


