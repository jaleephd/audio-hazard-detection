function decrease = calcSpectralDecrease(yn)

%function decrease = calcSpectralDecrease(y)
% yn (time,freq) is the normalised amplitudes of the frequency component (0-1)
%% y (time,freq) is the frequency component, magnitude and phase
%
% spectral decrease represents the amount of decrease in the spectral
% amplitude over the frequency spectrum. This is defined as:
%   decrease = (1 / sum_2:k ampl(k)) x sum_2:k (a(k)-a(1)) / (k-1)
% where k is the frequency components, a(1) is the 0 Hz (DC component)
%
% note this can be done with normalised (Yn) or unnormalised amplitudes (Y)
%      as the result is normalised anyway (just remove the abs() for yn)


decrease=zeros(size(yn,1),1); % for each time slice
k=size(yn,2); % the number of frequency components

for (i=1:size(yn,1))
  ampl=yn(i,:);
  a=sum((ampl(2:k)-ampl(1))/(k-1));
  total=sum(ampl(2:k));
  decrease(i)=a/total;
end


%decrease=zeros(size(y,1),1); % for each time slice
%k=size(y,2)); % the number of frequency components
%
%for (i=1:size(y,1))
%  ampl=abs(y(i,:));
%  a=sum((ampl(2:k)-ampl(1))/(k-1));
%  total=sum(ampl(2:k));
%  decrease(i)=a/total;
%end

