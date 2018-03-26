function [cnt,idx,width] = findPeaks(x,threshold)
%
% find the peaks in the signal, where a peak is a contiguous
% set of samples exceeding the threshold
%
% params:
%   - x contains a vector of signal samples
% returns:
%   - cnt is the number of peaks in the signal
%   - idx is the start of each peak
%   - width is the width of each peak


cnt=0;
xidx=zeros(size(x));
xwidth=zeros(size(x));


i=1;
while i<length(x)
  if x(i)>threshold
    cnt=cnt+1; % found the start of a peak
    xidx(cnt)=i;
    % skip to end of peak
    for j=i+1:length(x)+1 % +1 to take care of end of array condition
      if j>length(x) || x(j)==0
          break;
      end
    end
    xwidth(cnt)=j-i+1;
    i=j+1; % skip to next
  else
    i=i+1;
  end
end

idx=xidx(1:cnt);
width=xwidth(1:cnt);


