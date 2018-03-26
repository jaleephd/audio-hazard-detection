function [corr] = autocorr(Y)
% autocorr takes a row vector containing the FFT transform of a window 
% of temporal data, and returns a row vector of the same length
% containing the positive (first 1/2) and negative autocorrelation coefficients
%
% autocorrelation is defined as the correlation function applied to the same
% signal. ie: r=Corr(g,g)
% the result, r, of the autocorrelation is a vector of the coefficients,
% with the index giving the lag (time difference), with lags 0-N/2 in the
% first half of the vector, and the negative lags stored in elements
% from N-1 down. Note that r(1) is the correlation at zero lag
% (due to matlab's indexing)
%
% see also matlab functions: xcorr and corrcoef


% Ref: from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING, pp 545-546
%      http://hebb.mit.edu/courses/9.29/2002/readings/c13-2.pdf
%
% autocorrelation can be performed using the FFT
% as correllation coeeficient j for g(t) and h(t) is one member of the FFT pair
%      corr(g,h)_j = G_k H_k*
%    where g and h are time domain signals, of length N
%          * indicates the complex conjugate,
%          k is the frequency component
%
% To obtain the autocorrelation coefficients using the FFT:
% - take the FFT of a window of N samples of a signal (x)
%   to give N frequency components, X, 1/2 of which are -ve
% - multiply X by X* (the complex conjugate of X)
% - take the inverse FFT of the result
% - the resulting vector, r, of length N contains complex numbers 
%   with all imaginary parts equal to zero (which can be dropped).
%   The components or r are the values of the autocorrelation at different
%   lags, with lags 0-N/2 in the first half of the vector, and the negative
%   lags stored in elements from N-1 down

  % multiply Y by the complex conjugate of Y (Y*)
  r=Y.*conj(Y);

  % take the inverse FFT and keep only the real parts
  corr=real(ifft(r));





% in the time domain, correlation of g(t) and h(t) is defined as:
%   corr(g,h)_j = sum_k=0,N-1 g_j+k h_k
% where g and h are time domain signals, of length N
%   and j is the index of the coefficient with lag j

% the following code snippet to implement this is taken from Matlab Central:
%     speech processing tool, MEKHMOUKH Abdenour
%     http://www.mathworks.com.au/matlabcentral/fileexchange/18497

if 0
  % x is a series of N samples of a signal (amplitudes in time)
  N = size(x,2);
  x1 = x;
  x1(1,N+1:2*N-1) = zeros(1,N-1);
  x2 = zeros(N-1, 2*N-1);

  for i = 1:N-1
    x2(i,i:i+N-1) = x;
  end

  corr = 1/N*x1*x2';

end

