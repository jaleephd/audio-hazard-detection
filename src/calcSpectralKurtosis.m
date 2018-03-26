function flatness = calcSpectralKurtosis(yn,freqs,centroids)

% yn (time,freq) is normalised amplitudes of the frequency component (0-1)
% freqs is a list of frequencies in the spectrum (or band)
% centroids is a list of the centroids for each time slice
%
% the spectral kurtosis, or flatness, is a measure of the flatness of a
% distribution around its mean value, it is defined as:
% kurtosis = sum ((k-Cf)^4 * |Xn(k)|) / Sf^4
%   where k is the (center?) frequency of each component,
%         Cf is the centroid
%         |Xn(k)| is the normalised frequency component magnitude
%         Sf is the spread, found using
%           Sf = sqrt(sum((k-Cf))^2 * |X(k)|)


% calculate center freq of each band
f=freqs;
fhi=-1;
for i=1:length(f)
  flo=fhi+1;
  fhi=freqs(i);
  f(i)=(fhi+flo)/2;
end

flatness=zeros(size(yn,1),1); % for each time slice

for (i=1:size(yn,1))
  a=sum(((f-centroids(i)).^4).*yn(i,:));
  spread=sqrt(sum(((f-centroids(i)).^2).*yn(i,:)));
  flatness(i)=a/spread^4;
end

