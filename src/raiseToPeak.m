function y = raiseToPeak(x)
%
% sets all values within a non-zero set of consecutive samples to the maximum
% (peak) value within that set
%
%
% params:
%   - x contains a vector of signal samples
% returns:
%   - y is the modified signal


y=zeros(size(x));

i=1;
while i<length(x)
  if x(i)>0
    peakx=0;
    % find end of consecutive set of non-zero samples
    for j=i+1:length(x)+1 % +1 to take care of end of array condition
      if j>length(x) || x(j)==0
          break;
      else
        if x(j)>peakx
          peakx=x(j);
        end
      end
    end
    y(i:j-1)=peakx;  % set all consecutive samples to the peak value
    i=j+1; % skip to next sample
  else
    i=i+1;
  end
end


