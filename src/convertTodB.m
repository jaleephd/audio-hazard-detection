function Yd = convertTodB(Y)

% get amplitude and convert to dB
% taken from
%   http://hans.fugal.net/blog/2009/04/02/dft-magnitudepower-spectra
%  and
%   myspecgram.m by Paul Kienzle; modified by Sean Fulop March 2002

Yd = 20*log10(abs(Y));

% uncomment for debugging/visualisation
% disp('click graph to continue')
% k = waitforbuttonpress(); % pause on keystroke or mouse click
% bar(fr,Y);

