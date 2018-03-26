function ste = shortTimeEnergy(frame)
% calculate short-time energy for a time frame
% by taking the average of the the squared samples in the frame
% (note can apply a Hamming windown first for smoothing)
% STE = 1/N sum_n=1_N x(n)^2

ste = 1/length(frame) * sum(abs(frame.^2));

