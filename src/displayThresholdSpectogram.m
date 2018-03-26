function h = displayThresholdSpectogram(yd, fs, step, winLength, threshold)

% display freq comonents if signal > threshold dB : Y axis Freq, X axis time  
% yd contains:
%    - each row represents a time slice
%    - each column represents a frequency component
%    yd (time,freq) = amplitude (or mean) of that frequency component in dB
% fs is the sampling frequency (Hz), or an array of frequencies
% step is how often STFT was performed (in seconds)
% winLength is the length of the SFTT window (in samples)
% threshold is in dB, and can either be a single value, a vector (one for each
%           time slice), or a matrix thresholds, for each time slice, freq band

if (isscalar(fs))
  nf=floor((winLength+1)/2); % number of freq components (Nyquist)
  fr = [0:nf-1]*fs/winLength;  % generate the frequency axis for plotting
else
  fr=fs;
end

% find freq components (at given times) that have an amplitude above the threshold
if size(threshold,1)==1
  [t,f]=find(yd>=threshold); 
elsif size(threshold,2)==1 % a vector of thresholds
  t=[];
  f=[];
  for i=1:size(yd,1) % for each time slice
    [ti,fi]=find(yd(i,:)>=threshold(i));
    t=[t ti*i];
    f=[f fi];
  end
else % a matrix of thresholds, one for each time and frequency component
  [t,f]=find(yd-threshold>=0); 
end


tv=(t-1)*step; % convert from time indexes to time in seconds
fv=fr(f); % convert from freq indexes to frequency in Hz


clf reset % clear graph

% setup graph properties (log scale on Y axis for plotting with semilogy)
h=gca;
if (isscalar(fs))
  set(h,'yscale','log');
  ylim([10 fr(end)]); % ignore freq of 0
  ylabel('Frequency (Hz)');
else
  ylabel('Frequency Band (index)');
end

xlim([0 floor(tv(end)+1)]);
grid on;

title('Audio Frequency Components (thresholded)');
xlabel('time (sec)');

hold on;

if (isscalar(fs))
  semilogy(tv(:),fv(:),'.');
  %scatter(tv(:),fv(:),'filled');
else
  scatter(tv(:),f(:),'filled');
  %plot(tv(:),f(:),'.');
end


