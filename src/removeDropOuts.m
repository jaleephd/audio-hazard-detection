function y = removeDropOuts(x, width)
%
% removes drop outs to 0 if width of drop out is less than the specified width
%
% params:
%   - x contains a vector of signal samples
%   - width determines the narrowest valley not considered a drop out
% returns:
%   - y is the cleaned signal, with drop outs filled with the average of
%       the surrounding values



y=x;
i=1;
while i<length(x)
  if x(i)==0
    for j=i+1:length(x)+1 % +1 to take care of end of array condition
      if j>length(x) || x(j)>0
          break;
      end
    end
    if j-i+1 < width % if it's a drop out
      %sprintf('i=%d j=%d - %f - %f sec',i,j,i*0.050,j*0.050)
      if i==1 % at start of signal samples
        if j>length(x)
          fill=0; % this shouldn't happen
        else
          fill=x(j); % fill with following sample value
        end
      else % not at start
        if j>length(x) % if at end
          fill=x(i-1); % fill with previous sample value
        else
          fill=(x(i-1)+x(j))/2; % fill with average of prev & next samples
        end
      end
      y(i:(min(length(x),j)))=fill;
    end
    i=j+1; % skip to next
  else
    i=i+1;
  end
end



