function zcr = zeroCrossingRate(frame)
% ZCR = 1/2N sum_i=2:N |sgn(x(i))-sgn(x(i-1)|
% where sgn(x(i)) is +1 if x(i) >= 0, -1 if x(i)<0

sgn = (frame>=0) + (-1)*(frame<0);
sgn2 = [sgn(2:end); 0]; % add dummy on the end so same length
zcr = 1/(2*length(frame)) * sum(abs(sgn2-sgn));


