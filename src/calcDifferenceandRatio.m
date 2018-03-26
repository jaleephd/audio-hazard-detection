function [diff,diffratio] = calcDifferenceandRatio(a)

% a is a vector of recorded feature values that vary over time
% the difference between consecutive samples is calculated and returned in diff
% the ratio of the difference / (curr + prev samples) is returned in diffratio
% Note that diff and diffratio are one element shorter than vector a

diff=a(2:end)-a(1:end-1);
diffratio=diff./(a(1:end-1)+a(2:end));

