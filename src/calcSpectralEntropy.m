function entropy = calcSpectralEntropy(yn,freqs)

% yn (time,freq) is normalised amplitudes of the frequency component (0-1)
% freqs is a list of frequencies in the spectrum (or band)

% http://www.dsprelated.com/showmessage/108326/1.php
% Spectral entropy is calculated as:
% 1) Normalize the power spectrum:
%    Q(f)=P(f)/sum(P(f))  (where P(f) is the power spectrum)
%
% 2) Transform with the Shannon function:
%    H(f)=Q(f)[log(1/Q(f))]
%
% 3) Spectral entropy:
%    E=sum(H(f))/log(Nf)  (where Nf is the number of frequency components.
%
% For wavelet spectral entropy, replace step (2) with:
% 2) Transform with the Shannon function: 
%    Hi=Pi[log(1/Pi)]
% 3) Wavelet spectral entropy:
%    E=sum(Hi)/log(Ni) (where Ni is the number of subbands)


% power spectrum is already normalised (yn)
% E=sum (yn*log(1/yn)/log(length(freqs))

lognf=log(length(freqs));
entropy=zeros(size(yn,1),1); % for each time slice
for (i=1:size(yn,1))
  Hf=yn(i,:).*-log(yn(i,:));
  e=sum(Hf)/lognf;
  entropy(i) = e;
end


% En(i) = -sum(x.*log2(x));

