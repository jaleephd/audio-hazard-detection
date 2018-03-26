function energy = calcSpectralEnergy(y)
% y is the FFT transformed signal (or a sub band of it), and contains
%   the magnitude and phase of spectral components (complex number),
%   indexed as:
%    - each row represents a time slice
%    - each column represents a frequency component

% the energy of a frequency band is calculated as the sum of magnitude^2
% of the frequency components in the band:
%    En,m = sum k=kfirst_klast |Xk|^2
%    where |Xk| is the magnitude of the kth component of the FFT of the
%    (nth) windowed frame (ie the time axis)
%    kfirst and klast are the first and last FFT components belonging to
%    the (mth) frequency band


energy=zeros(size(y,1),1); % for each time slice
for i=1:size(y,1)
  energy(i)=sum(abs(y(i,:)).^2); % sum for all frequency components
end


