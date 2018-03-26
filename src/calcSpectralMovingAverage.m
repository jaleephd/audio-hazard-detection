function yma = calcSpectralMovingAverage(yd, n)

% yd contains:
%    - each row represents a time slice
%    - each column represents a frequency component
%    yd (time,freq) = amplitude of that frequency component in dB
% n is length of moving average

yma=zeros(size(yd));

% this uses the matlab function y=filter(b,a,x)
% where 
%      a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                            - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%   If a(1) is not equal to 1, FILTER normalizes the filter
%   coefficients by a(1). 
% eg for length 4 window: a=1, b = [1/4 1/4 1/4 1/4]
%    y(n)=1/4 x(n) + 1/4 x(n-1) + 1/4 x(n-2) + 1/4 x(n-3)

a=1;
b=ones(1,n) * 1/n;

% this could be optimised for matrix, but for clarity, do this way
for i=1:size(yd,2) % for each frequency component (column)
  yma(:,i)=filter(b,a,yd(:,i)); % calculate moving average
end
  

