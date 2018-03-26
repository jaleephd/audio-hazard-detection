function ZCR=calcZeroCrossingRate(x,winStep,winLength)

% create Hamming window
hWin = hamming(winLength);

nWindows=getNumWindows(x,winStep,winLength);

ZCR = zeros(nWindows,1);
for i=1:nWindows
    xoffset = ((i-1)*winStep)+1;
    frame = x(xoffset:xoffset+winLength-1).*hWin; % apply a Hamming window to this time frame
    ZCR(i) = zeroCrossingRate(frame);
end




