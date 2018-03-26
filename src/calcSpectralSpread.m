function spread = calcSpectralSpread(yn,freqs,centroids)

% yn (time,freq) is normalised amplitudes of the frequency component (0-1)
% freqs is a list of frequencies in the spectrum (or band)
% centroids is a list of the centroids for each time slice
%
% the spectral spread (bandwidth) is defined as the spread of the
% spectrum around its mean value (ie the variance).
% this is calculated as the sqrt of the average of the difference
% between the (center?) freq of each component and the spectral centroid
% ie: spread = sqrt(sum((f(n)-sc))^2 y(n))
%     where y(n) represents the normalised magnitude of bin number n,
%           f(n) represents the center frequency of that bin.
%           sc represents the spectral centroid


% calculate center freq of each band
f=freqs;
fhi=-1;
for i=1:length(f)
  flo=fhi+1;
  fhi=freqs(i);
  f(i)=(fhi+flo)/2;
end

spread=zeros(size(yn,1),1); % for each time slice
for (i=1:size(yn,1))
  spread(i)=sqrt(sum(((f-centroids(i)).^2).*yn(i,:)));
end

