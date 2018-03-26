function y = removeSpikes(x, width, thresh)
%
% removes spikes emanating from 0 and reaching or exceeding thresh
%         if spikes are less than the specified width
%
% params:
%   - x contains a vector of signal samples
%   - width determines the narrowest peak not considered a spike
%   - thresh is the minimum amplitude spike from 0 to be affected
%     and determines the value of the despiked signal (spikes are
%     thresholded to this value)
% returns:
%   - y is the cleaned signal, with spikes thresholded to the given thresh

% calculate width of activity peaks/spikes

xwidth=zeros(size(x));
i=1;
while i<length(x)
  if x(i)>0
    for j=i+1:length(x)+1 % +1 to take care of end of array condition
      if j>length(x) || x(j)==0
          break;
      end
    end
    xwidth(i:j-1)=j-i+1;
    i=j+1; % skip to next
  else
    i=i+1;
  end
end


% eliminate spikes in activity caused by noise
% that is spikes eminating from 0 and at least reaching the threshold
% and then dropping back to zero
y=x;
for i=1:length(y)
  if xwidth(i)<width
    y(i)=min(thresh,x(i));
  end
end


