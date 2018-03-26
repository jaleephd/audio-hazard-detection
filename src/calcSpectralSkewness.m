function skew = calcSpectralSkew(yn,freqs,centroids)

% yn (time,freq) is normalised amplitudes of the frequency component (0-1)
% freqs is a list of frequencies in the spectrum (or band)
% centroids is a list of the centroids for each time slice
%
% Spectral skewness describes the asymmetry of the frequency distribution
% and is is defined as:
%   skew = sum (k-Cf)^3 * |Xn(k)| / Sf^3
%   where k is the (center?) frequency of each component,
%         Cf is the centroid
%         |Xn(k)| is the normalised frequency component magnitude
%         Sf is the spectral spread (variance about mean), found using
%            Sf = sqrt(sum((k-Cf))^2 * |Xn(k)|)


% calculate center freq of each band
f=freqs;
fhi=-1;
for i=1:length(f)
  flo=fhi+1;
  fhi=freqs(i);
  f(i)=(fhi+flo)/2;
end

skew=zeros(size(yn,1),1); % for each time slice

for (i=1:size(yn,1))
  a=sum(((f-centroids(i)).^3).*yn(i,:));
  spread=sqrt(sum(((f-centroids(i)).^2).*yn(i,:)));
  skew(i)=a/spread^3;
end

