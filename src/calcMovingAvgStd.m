function [m,s] = calcMovingAvgStd(x, n)
%
% calculate the moving average and standard deviation of signal x
% params:
%   - x contains a vector of signal samples
%   - n is length of moving average window
% returns:
%   - m the moving average (a vector the same length as x)
%   - s the moving standard deviation (vector same length as x)


% calculate moving average using the matlab function y=filter(b,a,x)
% where 
%      a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                            - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%   If a(1) is not equal to 1, FILTER normalizes the filter
%   coefficients by a(1). 
% eg for length 4 window: a=1, b = [1/4 1/4 1/4 1/4]
%    y(n)=1/4 x(n) + 1/4 x(n-1) + 1/4 x(n-2) + 1/4 x(n-3)

a=1;

b=ones(1,n);
m=filter(b,a,x) * (1/n); % calculate moving average
 
% calculate moving standard deviation from prev n elements
% using std = sqrt((sum(x.^2) - n*xbar.^2)/(n-1))
% ref: http://www.mathworks.com/matlabcentral/fileexchange/9428
%s=sqrt((filter(b,a,x2) - (filter(b,a,x).^2)*(1/n))/(n-1));


% ref: http://newsgroups.derkeiler.com/Archive/Comp/comp.soft-sys.matlab/2006-05/msg00833.html
%  std = sqrt(1/(n-1)*sum((x-xm).^2))
%       = sqrt(n/(n-1)*(x2m - xm^2))
% where xm is (moving) average and x2m the (moving) average of its squares
x2=x.^2;
x2m=filter(b,a,x2)/n;
xm=filter(b,a,x)/n;
s=sqrt(n/(n-1)*(x2m-xm.^2));

% fix the first n elements
for i = 1:(n-1)
  s(i) = std(x(1:i));
end
%s(1:n-1)=std(x(1:n-1));
%s(1:n-1)=s(n);
%s(1:n-1)=0;


