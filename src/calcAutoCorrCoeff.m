function coeff = calcAutoCorrCoeff(y)
% y is the FFT transformed signal (or a sub band of it), and contains
%   the magnitude and phase of spectral components (complex number),
%   indexed as:
%    - each row represents a time slice
%    - each column represents a frequency component
%
% the auto correlation coefficients are calculated for each windowed
% sample, and the result being a row vector of the coefficients,
% with the index giving the lag (time difference), with lags 0-N/2 in the
% first half of the vector, and the negative lags stored in elements
% from N-1 down

coeff=zeros(size(y)); % for each time slice and coefficient
for i=1:size(y,1)
  coeff(i,:)=autocorr(y(i,:));
end



