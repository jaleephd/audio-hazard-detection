function flux = calcSpectralFlux(yn)

% yn (time,freq) is normalised amplitudes of the frequency component (0-1)
%
% spectral flux (spectral variation) is calculated as the normalised
% cross-correlation between two successive amplitude spectra
% flux = 1 - (sum f_t-1(n) f_t(n)) / (sqrt(sum f_t-1(n)^2) sqrt(sum f_t(n)^2))
%   where f_t(n) is the (normalised) amplitude of frequency component n
%                at time t

flux=zeros(size(yn,1),1); % for each time slice
flux(1)=0;
for (i=2:size(yn,1))
  a=sum(yn(i-1,:).*yn(i,:));
  b=sqrt(sum(yn(i-1,:).^2));
  c=sqrt(sum(yn(i,:).^2));
  flux(i)=1-(a*b/c);
end

