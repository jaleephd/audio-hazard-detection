function Y = getFreqComponents(XX, fs, n)

% XX is a vector (of length n) of (windowed) time domain data,
%    with a sampling freq fs
% Y is a vector of length n/2 of the signal's freq component amplitudes in dB

% perform DFT on window

Y = fft(XX); % contains + and -ve freq

%Y = fftshift(Y); % zero center freq
% extract +ve freq components
if (mod(n,2)==1)
    k=(n+1)/2;
else
    k=n/2;
end
Y=Y(1:k,:); % keep the +ve components

% uncomment for debugging/visualisation
% fr = [0:k-1]*fs/n;  % generate the frequency axis
% stem(fr,abs(Y));  %use abs command to get the magnitude
% %similary, use angle command to get the phase plot!

% don't do this yet, as will normalise amplitudes beforehand
%% take amplitude component and convert to dB
%Y=convertTodB(Y);


