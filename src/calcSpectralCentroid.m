function centroids = calcSpectralCentroid(yn,freqs)

% yn (time,freq) is normalised amplitudes of the frequency component (0-1)
% freqs is a list of frequencies in the spectrum (or band)
%
% spectral centroid is calculated as the weighted mean of the frequencies
% present in the signal, with their magnitudes as the weights
% C = (sum f(n) y(n)) / sum y(n)
%     where y(n) represents the (normalised) magnitude of bin number n,
%       and f(n) represents the center frequency of that bin.


% calculate center freq of each band
f=freqs;
fhi=-1;
for i=1:length(f)
  flo=fhi+1;
  fhi=freqs(i);
  f(i)=(fhi+flo)/2;
end

% calculate centroids
centroids=zeros(size(yn,1),1); % centroid for each time slice
for (i=1:size(yn,1))
  total=sum(yn(i,:));
  wmag=sum(f.*yn(i,:));
  centroids(i) = wmag/total;
end

