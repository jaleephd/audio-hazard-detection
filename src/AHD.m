%function FF = AHD(fileName, win, step)

clear;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% this section contains parameters that are changed from run to run for experiments


% the scenarios for testing, uncomment the one to be tested
%fileName='..//data//01_Beside ihbi intersection 2009-06-18 14-51-29.wav'
%fileName='..//data//03_Bus stop x busy road down suburban st 2009-06-19 17-53-55.wav'
%fileName='..//data//06_Cafe morning 2009-06-19 9-00-09.wav'
%fileName='..//data//08_Intersection - creek & adelaide 2009-06-19 9-36-52.wav'
%fileName='..//data//09_King George to Myer fri eve 2009-06-19 17-18-14.wav'
%fileName='..//data//11_Park 2009-06-18 15-15-14.wav'
%fileName='..//data//13_Suburban to main in rain 2009-06-22 8-12-05.wav'
%fileName='..//data//14_Woolies 2009-06-18 14-55-33.wav'
%fileName='..//data//38_Malbon st 1 2009-08-11 9-46-43.wav'
fileName='..//data//39_Malbon st 2 2009-08-11 9-49-01.wav'


% set this if want to save graph images using export_fig
saveimg=0 % save nothing
%saveimg=1 % save all graphs and results
%saveimg=-1 % only save method variant results, not spectral graphs


%IMPORTANT: sub band thresholding and scaling was done with moving average,
%      and SD of the moving average (over all time) as testing over limited time frames (80 sec)
%      this improves accuracy for identifying useful features, but isn't what would be done in a real implementation
%      if do over entire recordings (or implementing in real time), use the moving SD
%      by giving MsdWin a value > 0 
MsdWin = 0 % SD over entire recording (default)
%MsdWin = 30 % using a moving standard deviation - so set moving std dev window for scaling and thresholds


% determines whether roll off is used to influence sensitivity (adjust threshold)
% if set > 0 then roll off is used
userolloff=1 % use roll off (default)
%userolloff=0 % don't use roll off


% optionally define start and end of section, in seconds (use -1 or 0 for ignore)
%start=0;
%stop=0;
start=0
stop=80


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


win = 0.100; % 100 msec window
step = 0.050; % do STFT every 50 ms


% note that for mavWin = 10 and sdthresh = 2.5 gives eqiuv to spectogram 16 at threshold of ~ -60dB
% note that for mavWin = 5 and sdthresh = 3 gives eqiuv to spectogram 16 at threshold of ~ -70dB
mavWin = 10.0; % moving average window for smoothing mean for adaptive threshold
mavWin2 = 2.5; % moving average window for calculating SNR
mavWin3 = 10.0; % moving average window for sub band features
mavWin4 = 1.0; % moving average window for smoothing roll off

mWinLen=floor(mavWin/step);
mWinLen2=floor(mavWin2/step);
mWinLen3=floor(mavWin3/step);
mWinLen4=floor(mavWin4/step);
mStdWinLen=floor(MsdWin/step);

sdthresh=2.5; % threshold defined in terms of amplitude above smoothed mean, this is affected by mavWin
bandsdthresh=1;%0.80; % threshold for frequency bands, defined in terms of amplitude above smoothed mean
thresholddB = -30; % dB (depreceated)


% width of median filter to apply to signals for detecting vehicle activity
% this is to eliminate (noise) spikes. note there is a trade off
% between eliminating noise spikes and dampening the detection ability
medfiltw=floor(0.250/step); % 250ms / STFT sample rate    
    
% min width of a thresholded band activity peak to be considered
% less than this width is considered to be a noise spike
minvactw=floor(0.500/step); % 0.5 sec / STFT sample rate

% as above but used for results presentation
minvactresw=floor(1.000/step); % 1.0 sec / STFT sample rate

% if using a moving SD, the system will be a lot more sensitive to noise,
% so the following are scaled to make it less sensitive
if MsdWin>0
    sdthresh=sdthresh*1.5; % this just affects the spectogram display
    bandsdthresh=bandsdthresh*1.5;
    medfiltw=medfiltw*2;
    minvactw=minvactw*1.5;
end

% define the freq bands, noting that there is a 10Hz granularity limit, and
% that for lofi recordings at 22050 Hz, 11025 Hz is the upper limit
% each band is defined as from previous band high limit + 1 Hz (or 0 Hz for
% first band) to the (upper limit) frequency entered in the vector
%           1  2  3   4   5   6   7   8   9   10  11   12   13   14   15       16
freqbands=[40 70 110 150 200 250 300 400 500 750 1000 1500 2000 3000 5000]; % 11025

% load wav file
[x, fs] = wavread(fileName);

% if only want a section of the file, then just keep the defined section
if start>0 || stop>0
    if stop<=start || stop*fs>length(x)
        stop=floor(length(x)/fs);
    end
    if start<0
        start=0;
    end
    x=x((start*fs)+1:stop*fs);
end


winLength=floor(win*fs);
winStep=floor(step*fs);
nf=floor((winLength+1)/2); % number of freq components (Nyquist)
%fr = [0:nf-1]*fs/winLength;  % generate the frequency axis for plotting etc
fr = [0:nf-1]/win;  % generate the frequency axis for plotting etc: freq will be 0, 1/size of win, .. fs/2


% determine indexes (into fr, and hence FFT) of frequency bands
j=1;
freqbandsidx=zeros(size(freqbands)); % by default bands contain nothing
for i=1:length(fr);
    if fr(i)>freqbands(j)
        freqbandsidx(j)=i-1; % set index to be that of highest freq within band
        j=j+1;
        if j>length(freqbands) % freqencies outside of what were defined (were interested in)
            freqbands=[freqbands fr(end)]; % append freq band that covers the remaining frequencies
            freqbandsidx(j)=length(fr); % last band ends at highest recorded frequency
            break;
        end
    end    
end


% if bands exceed recorded frequencies, remove unused bands
if j<length(freqbands)
    freqbands=freqbands(1:j);
end

% create vector of freq band sizes in terms of number of frequency indexes
% each contains, for creating cell array of bands of spectral data
freqbandsizes=zeros(size(freqbandsidx));
hiidx=0;
for i=1:length(freqbandsidx)
    loidx=hiidx+1;
    hiidx=freqbandsidx(i);
    freqbandsizes(i)=hiidx-loidx+1;
end


nWindows=getNumWindows(x,winStep,winLength);
y=zeros(nWindows,nf); % this will contain the amplitudes of the freq components for the entire recording

% time vectors for plotting graphs
ts=(0:length(x)-1) / fs; % original (based on sampling freq)
tsw=(0:nWindows-1) * step; % for windowed measures



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- perform STFT --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('performing FFT ...')

% perform STFT
% y contains magnitude and phase of spectral components (complex number),
%   and is indexed as:
%    - each row represents a time slice
%    - each column represents a frequency component
y=calcSTFT(x,winStep,winLength);


% normalise signal ao that maximum amplitude will be 0 dB, and convert
% amplitudes to dB, clipping signals with amplitude less than -100 dB

% yd (time,freq) = amplitude of that frequency component in dB
% yn (time,freq) = normalised amplitude of that frequency component (0-1)
[yd,yn,maxampl,minampl]=normaliseTodB(y,0,-1,-100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- calculate temporal features --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate short-time energy for a time frame (i)
% by taking the average of the the squared samples in the frame
% (note we apply a Hamming windown first for smoothing)
STE = calcShortTimeEnergy(x,winStep,winLength);

% calculate zero-crossing rate
ZCR = calcZeroCrossingRate(x,winStep,winLength);

% calculate RMS (ie power)
rms=calcRMS(x,winStep,winLength);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- calculate spectral features --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate auto correlation coefficients
%    - each row represents a windowed time slice
%    - each column represents a correlation coefficient
%      0-n/2 are +ve (incl 0), n-1 - n/2+1 are -ve
acoeffs=calcAutoCorrCoeff(y);

% calculate spectral/signal energy
energy=calcSpectralEnergy(y);

% calculate centroid
centroid=calcSpectralCentroid(yn,fr);

% calculate bandwidth
bandwidth=calcSpectralSpread(yn,fr,centroid);

% calculate entropy
entropy=calcSpectralEntropy(yn,fr);

% calculate spectral flux
flux=calcSpectralFlux(yn);

% calculate spectral roll off
[rolloff,rolloffidx]=calcSpectralRollOff(yn,fr);

% calculate spectral flatness (kurtosis)
flatness=calcSpectralKurtosis(yn,fr,centroid);

% calculate spectral skewness
skewness=calcSpectralSkewness(yn,fr,centroid);

% calculate spectral decrease
spdecrease=calcSpectralDecrease(yn);
%spdecrease=calcSpectralDecrease(y);


% smooth signal with a median filter (1 sec window)
sSTE = medfilt1(STE,mWinLen4);
sZCR = medfilt1(ZCR,mWinLen4);
srms = medfilt1(rms,mWinLen4);
sspdecrease = medfilt1(spdecrease,mWinLen4);
sskewness = medfilt1(skewness,mWinLen4);
sflatness = medfilt1(flatness,mWinLen4);
sflux = medfilt1(flux,mWinLen4);
scentroid = medfilt1(centroid,mWinLen4);
srolloff = medfilt1(rolloff,mWinLen4);
sentropy = medfilt1(entropy,mWinLen4);
sbandwidth = medfilt1(bandwidth,mWinLen4);
senergy = medfilt1(energy,mWinLen4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract some useful features from roll off
%srolloff = medfilt1(rolloff,mWinLen4); % smoothed rolloff (1 sec window)
[rolloffMA,md] = calcMovingAvgStd(srolloff,mWinLen4);

% 2 different ways of calculating slope
roff=(srolloff-rolloffMA)./srolloff;
slopeRO=gradient(srolloff);

% calculate the logarithmic amount of change in the last mWinLen second window
roffD=zeros(size(rolloff));
for i=1:length(rolloff)
  j=max(1,i-mWinLen+1);
  roffD(i)=log10(max(srolloff(j:i)))-log10(min(srolloff(j:i)));
end


if userolloff>0
  % use logarithmic (log10) changes in roll off to increase or decrease sensitivity 
  % to sub band amplitues and features (in following tests)
  %
  % roffD values typically range from <0.1 to >1.0
  % with higher values generally associated with vehicle passing
  % try simple approach:
  %   - values >= 0.65 increase likelyhood that is a vehicle (set threshold lower)
  %   - values 0.35 - 0.65 have no effect (leave threshold as is)
  %   - values < 0.35 reduce likelyhood that is a vehicle (set threshold higher)

  % give factor of 0.75, 1, 1.5 based on value of rofD
  % this will be used to scale threshold of band amplitude and features
  rofactor=((roffD<0.35)*0.5)+((roffD<0.65)*0.25)+0.75;

else
  rofactor=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%% graph temporal and global spectral features %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%if 0

clf reset % clear graph
hold on
title('x, STE, ZCR, flux');
xlabel('time (sec)');
ylabel('Scaled Amplitude');
ylim([-1 1]);

plot(ts(:),x(:),'g:');
plot(tsw(:),STE(:)/max(STE),'r-');
plot(tsw(:),ZCR(:)/max(ZCR),'b-');
plot(tsw(:),flux(:)/(max(flux)-min(flux)),'m-');

lgn=legend('amplitude','STE','ZCR', 'flux', 3);
h=gca;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
if saveimg>0
  export_fig('global_x+ste+zcr+flux.jpg','-jpg',h)
end


clf reset % clear graph
hold on
title('x, flatness, skewness, specdec, energy');
xlabel('time (sec)');
ylabel('Scaled Amplitude');
ylim([-1 1]);

plot(ts(:),x(:)/max(x),'g:');
plot(tsw(:),flatness(:)/max(flatness),'b-');
plot(tsw(:),skewness(:)/max(skewness),'c-'); % this follows flatness
plot(tsw(:),spdecrease(:)/(max(spdecrease)-min(spdecrease)),'m-');
plot(tsw(:),energy(:)/max(energy),'r-'); % this also follows signal amplitude
%plot(tsw(:),rms(:)/max(rms),'c-'); % this just follows signal amplitude

lgn=legend('amplitude','flatness', 'skewness', 'decrease', 'energy', 3);
h=gca;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
if saveimg>0
  export_fig('global_x+flat+skew+spdec+energy.jpg','-jpg',h)
end


clf reset % clear graph
hold on
title('x, bandwidth, centroid, entropy');
xlabel('time (sec)');
ylabel('Scaled Amplitude');
ylim([-1 1]);

plot(ts(:),x(:)/max(x),'g:');
plot(tsw(:),bandwidth(:)/max(bandwidth),'r-'); % this just follows signal amplitude
plot(tsw(:),centroid(:)/max(centroid),'b-');
plot(tsw(:),entropy(:)/max(entropy),'m-');

lgn=legend('amplitude','bandwidth', 'centroid', 'entropy', 3);
h=gca;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
if saveimg>0
  export_fig('global_x+bw+centroid+entropy.jpg','-jpg',h)
end


%end


%%%%%%%%%%%%%%%%%%%%%%%%%% graph smoothed temporal and global spectral features %%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 0
clf reset % clear graph
hold on
title('smoothed energy, BW, entropy, flatness, flux, spdec');
xlabel('time (sec)');
ylabel('Scaled Magnitude');
grid on;

plot(tsw(:),senergy(:)/max(senergy),'k-');
%plot(tsw(:),sSTE(:)/max(sSTE),'k-'); % exactly the same as energy
%plot(tsw(:),srms(:)/max(srms),'k-'); % magnified version of energy
plot(tsw(:),sbandwidth(:)/max(sbandwidth),'r-'); % this just follows signal amplitude
plot(tsw(:),sentropy(:)/max(sentropy),'b-');
plot(tsw(:),sflatness(:)/max(sflatness),'g-');
plot(tsw(:),sflux(:)/(max(flux)-min(flux)),'m-');
plot(tsw(:),sspdecrease(:)/(max(sspdecrease)-min(sspdecrease)),'c-');

lgn=legend('energy','BW', 'entropy', 'flatness', 'flux', 'spdec', 4);
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
end


%%%%%%%%%%%%%%%%%%%%%%%%%% graph velocity & accel of global signal amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%


% look at velocity and acceleration of signal amplitude (using RMS)
[drms, drrms] = calcDifferenceandRatio(rms);
[ddrms, ddrrms] = calcDifferenceandRatio(drms);

if 0
clf reset % clear graph
hold on
title('Signal Power, Velocity and Acceleration');
%title('Normalised Signal Power, Velocity and Acceleration');
xlabel('time (sec)');
ylabel('Amplitude');

plot(tsw(:),rms(:),'k-');
plot(tsw(:),[0; drms],'bo:');
plot(tsw(:),[0; 0; ddrms],'rx:');

%plot(tsw(:),rms(:)/max(rms),'k-');
%plot(tsw(:),[0; drms]/max(abs(drms)),'bo:');
%plot(tsw(:),[0; 0; ddrms]/max(abs(ddrms)),'rx:');

lgn=legend('mean','vel','acc',4);
%grid on;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graph Autocorrelation Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if 0
clf reset % clear graph
hold on
title('Autocorrelation Coefficients');
xlabel('time (sec)');
ylabel('Amplitude');
%ylim([-1 1]);
grid on;

plot(tsw(:),rms(:)/max(rms),'k:'); % this represents signal amplitude
%plot(ts(:),x(:)/max(x),'g:');

% each coeff scaled
%plot(tsw(:),acoeffs(:,1)/(max(acoeffs(:,1))-min(acoeffs(:,1))),'b:');
%plot(tsw(:),acoeffs(:,64)/(max(acoeffs(:,64))-min(acoeffs(:,64))),'r:');
%plot(tsw(:),acoeffs(:,128)/(max(acoeffs(:,128))-min(acoeffs(:,128))),'g:');
%plot(tsw(:),acoeffs(:,256)/(max(acoeffs(:,256))-min(acoeffs(:,256))),'m:');
%plot(tsw(:),acoeffs(:,512)/(max(acoeffs(:,512))-min(acoeffs(:,512))),'c:');

% scaled by 0th coefficient
plot(tsw(:),acoeffs(:,1)/(max(acoeffs(:,1))-min(acoeffs(:,1))),'b-');
%plot(tsw(:),acoeffs(:,2)/(max(acoeffs(:,1))-min(acoeffs(:,1))),'r:');
plot(tsw(:),acoeffs(:,64)/(max(acoeffs(:,1))-min(acoeffs(:,1))),'r-');
%plot(tsw(:),acoeffs(:,4)/(max(acoeffs(:,1))-min(acoeffs(:,1))),'g:');
plot(tsw(:),acoeffs(:,128)/(max(acoeffs(:,1))-min(acoeffs(:,1))),'g-');
%plot(tsw(:),acoeffs(:,8)/(max(acoeffs(:,1))-min(acoeffs(:,1))),'m:');
plot(tsw(:),acoeffs(:,256)/(max(acoeffs(:,1))-min(acoeffs(:,1))),'m-');
%plot(tsw(:),acoeffs(:,16)/(max(acoeffs(:,1))-min(acoeffs(:,1))),'c:');
plot(tsw(:),acoeffs(:,512)/(max(acoeffs(:,1))-min(acoeffs(:,1))),'c-');

% unscaled
%plot(tsw(:),acoeffs(:,1),'b:');
%plot(tsw(:),acoeffs(:,64),'r:');
%plot(tsw(:),acoeffs(:,128),'g:');
%plot(tsw(:),acoeffs(:,256),'m:');
%plot(tsw(:),acoeffs(:,512),'c:');

%lgn=legend('rms', 'coeff 0', 'coeff 1','coeff 3', 'coeff 7', 'coeff 15', 4);
lgn=legend('rms', 'coeff 0', 'coeff 63','coeff 127', 'coeff 255', 'coeff 511', 3);
h=gca;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
if saveimg>0
  export_fig('acoeffs.jpg','-jpg',h)
end

%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graph Spectral Roll Off %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if 0

% calculate the amount of change in the last mWinLen second window for graphing only
% this differs from roffD used in the heuristic
graffroffD=zeros(size(rolloff));
for i=1:length(rolloff)
  j=max(1,i-mWinLen+1);
  graffroffD(i)=max(srolloff(j:i))-min(srolloff(j:i));
end

% setup graph properties (log scale on Y axis for plotting with semilogy)
clf reset % clear graph
h=gca;
set(h,'yscale','log');
ylim([10 fr(end)]); % ignore freq of 0
%xlim([0 floor(ts(end)+1)]);
grid on;

title('Roll-off Frequency');
xlabel('time (sec)');
ylabel('Frequency (Hz)');

hold on;
semilogy(tsw(:),rolloff(:),'g--');
semilogy(tsw(:),rolloffMA(:),'b-');
semilogy(tsw(:),graffroffD(:),'r.-');

lgn=legend('rolloff','smoothed','change', 3);
h=gca;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
if saveimg>0
  export_fig('rolloff.jpg','-jpg',h)
end
%end



if 0
% setup graph properties (log scale on Y axis for plotting with semilogy)
clf reset % clear graph
h=gca;
set(h,'yscale','log');
ylim([10 fr(end)]); % ignore freq of 0
%xlim([0 floor(ts(end)+1)]);
grid on;

title('Roll-off Frequency');
xlabel('time (sec)');
ylabel('Frequency (Hz)');

hold on;
semilogy(tsw(:),rolloff(:),'b:');
%srolloff = medfilt1(rolloff,medfiltw);
srolloff = medfilt1(rolloff,mWinLen4);
semilogy(tsw(:),srolloff(:),'r-');
semilogy(tsw(:),rolloffMA(:),'g.-');
semilogy(tsw(:),abs(srolloff-rolloffMA),'m-');

lgn=legend('rolloff','smoothed','moving avg','deviation', 4);
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
end


if 0
% setup graph properties (log scale on Y axis for plotting with semilogy)
clf reset % clear graph
h=gca;
set(h,'yscale','log');
ylim([10 fr(end)]); % ignore freq of 0
%xlim([0 floor(ts(end)+1)]);
grid on;

title('Smoothed Roll-off Frequency Deviation');
xlabel('time (sec)');
ylabel('Frequency (Hz)');

hold on;
%semilogy(tsw(:),rolloff(:),'b:');
%srolloff = medfilt1(rolloff,medfiltw);
semilogy(tsw(:),srolloff(:),'r-');
%semilogy(tsw(:),rolloffMA(:),'g.-');

semilogy(tsw(:),roffD(:)*100,'m-');

lgn=legend('smoothed rolloff','windowed log deviation x 100', 4);
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
end


if 0
clf reset % clear graph
hold on
title('Roll-off Frequency derived Data');
xlabel('time (sec)');
ylabel('magnitude');

plot(tsw(:),roff(:)/(max(roff)-min(roff)),'b-');
plot(tsw(:),slopeRO(:)/max(slopeRO),'r-');

lgn=legend('smoothed-MA/smoothed','Gradient', 4);
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
end


if 0
clf reset % clear graph
hold on
grid on
title('Roll-off Frequency derived Data');
xlabel('time (sec)');
ylabel('magnitude');
plot(tsw(:),roffD(:),'m-');
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% display signal to noise ratios for global signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0

% calculate signal to noise ratio
snr_acoeff = calcSNR(acoeffs(:,1),mWinLen2);
snr_STE = calcSNR(STE,mWinLen2);
snr_ZCR = calcSNR(ZCR,mWinLen2);
snr_rms = calcSNR(rms,mWinLen2);
snr_spdecrease = calcSNR(spdecrease,mWinLen2);
snr_skewness = calcSNR(skewness,mWinLen2);
snr_flatness = calcSNR(flatness,mWinLen2);
snr_flux = calcSNR(flux,mWinLen2);
snr_centroid = calcSNR(centroid,mWinLen2);
snr_rolloff = calcSNR(rolloff,mWinLen2);
snr_entropy = calcSNR(entropy,mWinLen2);
snr_bandwidth = calcSNR(bandwidth,mWinLen2);
snr_energy = calcSNR(energy,mWinLen2);


clf reset % clear graph

% set up for logarithmic plotting on Y axis
h=gca;
set(h,'yscale','log');

hold on
title('SNR: STE, ZCR, flux, entropy');
xlabel('time (sec)');
ylabel('Signal to Noise Ratio');

semilogy(tsw(:),min(snr_STE(:),1000),'g-');
semilogy(tsw(:),min(snr_ZCR(:),1000),'r-');
semilogy(tsw(:),min(abs(snr_flux(:)),1000),'b-');
semilogy(tsw(:),min(snr_entropy(:),1000),'m-');

lgn=legend('STE','ZCR', 'flux', 'entropy',4);
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


clf reset % clear graph

% set up for logarithmic plotting on Y axis
h=gca;
set(h,'yscale','log');

hold on
title('SNR: roll off, centroid, flatness, skewness, specdec');
xlabel('time (sec)');
ylabel('Signal to Noise Ratio');

semilogy(tsw(:),min(snr_rolloff(:),1000),'g-');
semilogy(tsw(:),min(snr_centroid(:),1000),'r-');
semilogy(tsw(:),min(snr_flatness(:),1000),'b-');
semilogy(tsw(:),min(snr_skewness(:),1000),'m-');
semilogy(tsw(:),min(abs(snr_spdecrease(:)),1000),'c-');

lgn=legend('roll off','centroid', 'flatness', 'skewness', 'decrease', 4);
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


clf reset % clear graph

% set up for logarithmic plotting on Y axis
h=gca;
set(h,'yscale','log');

hold on
title('SNR: acoeff, rms, energy, bandwidth');
xlabel('time (sec)');
ylabel('Signal to Noise Ratio');

semilogy(tsw(:),min(snr_acoeff(:),1000),'g-');
semilogy(tsw(:),min(snr_rms(:),1000),'r-');
semilogy(tsw(:),min(snr_energy(:),1000),'b-');
semilogy(tsw(:),min(snr_bandwidth(:),1000),'m-');

lgn=legend('acoeff', 'rms', 'energy', 'bandwidth', 4);
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% graph Spectogram of all Frequency components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
% display freq comonents if signal > threshold dB : Y axis Freq, X axis time
h=displayThresholdSpectogram(yd,fs,step,winLength,thresholddB);
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
end

% generate a moving average of the magnitudes of the frequency components
% and use this average as the threshold, instead of a fixed threshold
% to determine the diffence between signal and moving average (ie noise)
% depreceated - it's basically a low pass filter
% NOTE: doesn't correlate well with what we are looking for
%yma=calcSpectralMovingAverage(yd,mWinLen); 
%ydd=abs(yd-yma);
%h=displayThresholdSpectogram(ydd,fs,step,winLength,10);
%h=displayThresholdSpectogram(yma,fs,step,winLength,-50);
%disp('click graph to continue')
%k = waitforbuttonpress(); % pause on keystroke or mouse click


% create a moving average and standard deviation of total signal over time,
% to look at spectral informaion that stands out, rather than just background noise
thresholds=calcSpectralThresholds(yd,mWinLen,mStdWinLen,sdthresh);


if 0
h=displayThresholdSpectogram(yd,fs,step,winLength,thresholds);
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% divide spectral data into frequency bands
% divide matrix into n cells of:
%   freq bands of various numbers of columns (frequencies),
%   each containing all the rows (ie over the entire time of the recording)
bands=mat2cell(y, size(y,1), freqbandsizes); % sub bands with freq components (phase and magnitude)
bandsdb=mat2cell(yd, size(yd,1), freqbandsizes); % sub bands with freq components in dB
bandsn=mat2cell(yn, size(yn,1), freqbandsizes); % sub bands with freq components in normalised magnitudes

% feature matrices

bandthresholds=zeros(size(yd,1), length(freqbands));
bandacoeffs=zeros(size(y,1), length(freqbands), 3); % keep 3 coeffs 
bandmean=zeros(size(yd,1), length(freqbands));
bandstd=zeros(size(yd,1), length(freqbands));
bandcentroid=zeros(size(yn,1), length(freqbands));
bandbandwidth=zeros(size(yn,1), length(freqbands));
bandentropy=zeros(size(yn,1), length(freqbands));
bandenergy=zeros(size(yn,1), length(freqbands));
bandenergyratio=zeros(size(yn,1), length(freqbands));
bandflatness=zeros(size(yn,1), length(freqbands));
bandskewness=zeros(size(yn,1), length(freqbands));
bandflux=zeros(size(yn,1), length(freqbands));
bandspectraldecrease=zeros(size(yn,1), length(freqbands));


% extracted features scaled in terms of standard deviation from moving average
scbandacoeffs=zeros(size(y,1), length(freqbands)); % only using 0th coeff 
scbandcentroid=zeros(size(yn,1), length(freqbands));
scbandbandwidth=zeros(size(yn,1), length(freqbands));
scbandentropy=zeros(size(yn,1), length(freqbands));
scbandenergy=zeros(size(yn,1), length(freqbands));
scbandenergyratio=zeros(size(yn,1), length(freqbands));
scbandflatness=zeros(size(yn,1), length(freqbands));
scbandskewness=zeros(size(yn,1), length(freqbands));
scbandflux=zeros(size(yn,1), length(freqbands));
scbandspectraldecrease=zeros(size(yn,1), length(freqbands));


% thresholded the extracted features based on standard deviation from moving average
thbandacoeffs=zeros(size(y,1), length(freqbands)); % only using 0th coeff 
thbandcentroid=zeros(size(yn,1), length(freqbands));
thbandbandwidth=zeros(size(yn,1), length(freqbands));
thbandentropy=zeros(size(yn,1), length(freqbands));
thbandenergy=zeros(size(yn,1), length(freqbands));
thbandenergyratio=zeros(size(yn,1), length(freqbands));
thbandflatness=zeros(size(yn,1), length(freqbands));
thbandskewness=zeros(size(yn,1), length(freqbands));
thbandflux=zeros(size(yn,1), length(freqbands));
thbandspectraldecrease=zeros(size(yn,1), length(freqbands));


% feature difference (velocity) matrices (note one element shorter in temporal dimension)
dbandacoeffs=zeros(size(y,1)-1, length(freqbands)); % only using 0th coeff 
dbandmean=zeros(size(yd,1)-1, length(freqbands));
dbandstd=zeros(size(yd,1)-1, length(freqbands));
dbandcentroid=zeros(size(yn,1)-1, length(freqbands));
dbandbandwidth=zeros(size(yn,1)-1, length(freqbands));
dbandentropy=zeros(size(yn,1)-1, length(freqbands));
dbandenergy=zeros(size(yn,1)-1, length(freqbands));
dbandenergyratio=zeros(size(yn,1)-1, length(freqbands));
dbandflatness=zeros(size(yn,1)-1, length(freqbands));
dbandskewness=zeros(size(yn,1)-1, length(freqbands));
dbandflux=zeros(size(yn,1)-1, length(freqbands));
dbandspectraldecrease=zeros(size(yn,1)-1, length(freqbands));

% feature difference (velocity) ratio matrices (note one element shorter in temporal dimension)
drbandacoeffs=zeros(size(y,1)-1, length(freqbands)); % only using 0th coeff 
drbandmean=zeros(size(yd,1)-1, length(freqbands));
drbandstd=zeros(size(yd,1)-1, length(freqbands));
drbandcentroid=zeros(size(yn,1)-1, length(freqbands));
drbandbandwidth=zeros(size(yn,1)-1, length(freqbands));
drbandentropy=zeros(size(yn,1)-1, length(freqbands));
drbandenergy=zeros(size(yn,1)-1, length(freqbands));
drbandenergyratio=zeros(size(yn,1)-1, length(freqbands));
drbandflatness=zeros(size(yn,1)-1, length(freqbands));
drbandskewness=zeros(size(yn,1)-1, length(freqbands));
drbandflux=zeros(size(yn,1)-1, length(freqbands));
drbandspectraldecrease=zeros(size(yn,1)-1, length(freqbands));

% feature acceleration matrices (note two elements shorter in temporal dimension)
ddbandacoeffs=zeros(size(y,1)-2, length(freqbands)); % only using 0th coeff 
ddbandmean=zeros(size(yd,1)-2, length(freqbands));
ddbandstd=zeros(size(yd,1)-2, length(freqbands));
ddbandcentroid=zeros(size(yn,1)-2, length(freqbands));
ddbandbandwidth=zeros(size(yn,1)-2, length(freqbands));
ddbandentropy=zeros(size(yn,1)-2, length(freqbands));
ddbandenergy=zeros(size(yn,1)-2, length(freqbands));
ddbandenergyratio=zeros(size(yn,1)-2, length(freqbands));
ddbandflatness=zeros(size(yn,1)-2, length(freqbands));
ddbandskewness=zeros(size(yn,1)-2, length(freqbands));
ddbandflux=zeros(size(yn,1)-2, length(freqbands));
ddbandspectraldecrease=zeros(size(yn,1)-2, length(freqbands));

% feature acceleration ratio matrices (note two elements shorter in temporal dimension)
ddrbandacoeffs=zeros(size(y,1)-2, length(freqbands)); % only using 0th coeff 
ddrbandmean=zeros(size(yd,1)-2, length(freqbands));
ddrbandstd=zeros(size(yd,1)-2, length(freqbands));
ddrbandcentroid=zeros(size(yn,1)-2, length(freqbands));
ddrbandbandwidth=zeros(size(yn,1)-2, length(freqbands));
ddrbandentropy=zeros(size(yn,1)-2, length(freqbands));
ddrbandenergy=zeros(size(yn,1)-2, length(freqbands));
ddrbandenergyratio=zeros(size(yn,1)-2, length(freqbands));
ddrbandflatness=zeros(size(yn,1)-2, length(freqbands));
ddrbandskewness=zeros(size(yn,1)-2, length(freqbands));
ddrbandflux=zeros(size(yn,1)-2, length(freqbands));
ddrbandspectraldecrease=zeros(size(yn,1)-2, length(freqbands));



%%%%%%%%%%%%%%%%%%%%%%%%%% compute spectral features for each band: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  
  loidx=freqbandsidx(i)-freqbandsizes(i)+1;
  hiidx=freqbandsidx(i);
  
  % calculate band features
  
  % calculate auto correlation coefficients
  %    - each row represents a windowed time slice
  %    - each column represents a correlation coefficient (there as many as
  %      freq components in the band). Note: all coefficients are positive
  %      as only kept +ve freq components
  bacoeffs=calcAutoCorrCoeff(bands{i});
  sz=size(bacoeffs,2); % number of coefficients

  bandacoeffs(:,i,1)=bacoeffs(:,1); % first coefficient
  bandacoeffs(:,i,2)=bacoeffs(:,floor(sz/2)+1); % middle coefficient
  bandacoeffs(:,i,3)=bacoeffs(:,end); % last coefficient

  bandthresholds(:,i)=calcSpectralThresholds(bandsdb{i},mWinLen,mStdWinLen,bandsdthresh*rofactor)'; % note scaling of threshold by rofactor
  %bandthresholds(:,i)=calcSpectralThresholds(bandsdb{i},mWinLen,mStdWinLen,bandsdthresh)';
  bandmean(:,i)=mean(bandsdb{i}')';
  bandstd(:,i)=std(bandsdb{i}')';
  bandcentroid(:,i)=calcSpectralCentroid(bandsn{i},fr(loidx:hiidx));
  bandbandwidth(:,i)=calcSpectralSpread(bandsn{i},fr(loidx:hiidx),bandcentroid(:,i));
  bandentropy(:,i)=calcSpectralEntropy(bandsn{i},fr(loidx:hiidx));
  bandenergy(:,i)=calcSpectralEnergy(bandsn{i});
  bandenergyratio(:,i)=energy./bandenergy(:,i);
  bandflatness(:,i)=calcSpectralKurtosis(bandsn{i},fr(loidx:hiidx),bandcentroid(:,i));
  bandskewness(:,i)=calcSpectralSkewness(bandsn{i},fr(loidx:hiidx),bandcentroid(:,i));
  bandflux(:,i)=calcSpectralFlux(bandsn{i});
  bandspectraldecrease(:,i)=calcSpectralDecrease(bandsn{i});

  % scale band features in terms of SD away from moving average
  scbandacoeffs(:,i)=scaleSignalToSD(bandacoeffs(:,i,1),mWinLen3,mStdWinLen);
  scbandcentroid(:,i)=scaleSignalToSD(bandcentroid(:,i),mWinLen3,mStdWinLen);
  scbandbandwidth(:,i)=scaleSignalToSD(bandbandwidth(:,i),mWinLen3,mStdWinLen);
  scbandentropy(:,i)=scaleSignalToSD(bandentropy(:,i),mWinLen3,mStdWinLen);
  scbandenergy(:,i)=scaleSignalToSD(bandenergy(:,i),mWinLen3,mStdWinLen);
  scbandenergyratio(:,i)=scaleSignalToSD(bandenergyratio(:,i),mWinLen3,mStdWinLen);
  scbandflatness(:,i)=scaleSignalToSD(bandflatness(:,i),mWinLen3,mStdWinLen);
  scbandskewness(:,i)=scaleSignalToSD(bandskewness(:,i),mWinLen3,mStdWinLen);
  scbandflux(:,i)=scaleSignalToSD(bandflux(:,i),mWinLen3,mStdWinLen);
  scbandspectraldecrease(:,i)=scaleSignalToSD(bandspectraldecrease(:,i),mWinLen3,mStdWinLen);

  % threshold band features based on SD away from moving average
  thbandacoeffs(:,i)=(abs(scbandacoeffs(:,i))>=bandsdthresh).*sign(scbandacoeffs(:,i));
  thbandcentroid(:,i)=(abs(scbandcentroid(:,i))>=bandsdthresh).*sign(scbandcentroid(:,i));
  thbandbandwidth(:,i)=(abs(scbandbandwidth(:,i))>=bandsdthresh).*sign(scbandbandwidth(:,i));
  thbandentropy(:,i)=(abs(scbandentropy(:,i))>=bandsdthresh).*sign(scbandentropy(:,i));
  thbandenergy(:,i)=(abs(scbandenergy(:,i))>=bandsdthresh).*sign(scbandenergy(:,i));
  thbandenergyratio(:,i)=(abs(scbandenergyratio(:,i))>=bandsdthresh).*sign(scbandenergyratio(:,i));
  thbandflatness(:,i)=(abs(scbandflatness(:,i))>=bandsdthresh).*sign(scbandflatness(:,i));
  thbandskewness(:,i)=(abs(scbandskewness(:,i))>=bandsdthresh).*sign(scbandskewness(:,i));
  thbandflux(:,i)=(abs(scbandflux(:,i))>=bandsdthresh).*sign(scbandflux(:,i));
  thbandspectraldecrease(:,i)=(abs(scbandspectraldecrease(:,i))>=bandsdthresh).*sign(scbandspectraldecrease(:,i));
  
  % calculate the change over time (velocity) of band features
  [dbandacoeffs(:,i), drbandacoeffs(:,i)] = calcDifferenceandRatio(bandacoeffs(:,i,1));
  [dbandmean(:,i), drbandmean(:,i)] = calcDifferenceandRatio(bandmean(:,i));
  [dbandstd(:,i), drbandstd(:,i) ] = calcDifferenceandRatio(bandstd(:,i));
  [dbandcentroid(:,i), drbandcentroid(:,i)] = calcDifferenceandRatio(bandcentroid(:,i));
  [dbandbandwidth(:,i), drbandbandwidth(:,i)] = calcDifferenceandRatio(bandbandwidth(:,i));
  [dbandentropy(:,i), drbandentropy(:,i)] = calcDifferenceandRatio(bandentropy(:,i));
  [dbandenergy(:,i), drbandenergy(:,i)] = calcDifferenceandRatio(bandenergy(:,i));
  [dbandenergyratio(:,i), drbandenergyratio(:,i)] = calcDifferenceandRatio(bandenergyratio(:,i));
  [dbandflatness(:,i), drbandflatness(:,i)] = calcDifferenceandRatio(bandflatness(:,i));
  [dbandskewness(:,i), drbandskewness(:,i)] = calcDifferenceandRatio(bandskewness(:,i));
  [dbandflux(:,i), drbandflux(:,i)] = calcDifferenceandRatio(bandflux(:,i));
  [dbandspectraldecrease(:,i), drbandspectraldecrease(:,i)] = calcDifferenceandRatio(bandspectraldecrease(:,i));

  % calculate the acceleration of band features
  [ddbandacoeffs(:,i), ddrbandacoeffs(:,i)] = calcDifferenceandRatio(dbandacoeffs(:,i));
  [ddbandmean(:,i), ddrbandmean(:,i)] = calcDifferenceandRatio(dbandmean(:,i));
  [ddbandstd(:,i), ddrbandstd(:,i) ] = calcDifferenceandRatio(dbandstd(:,i));
  [ddbandcentroid(:,i), ddrbandcentroid(:,i)] = calcDifferenceandRatio(dbandcentroid(:,i));
  [ddbandbandwidth(:,i), ddrbandbandwidth(:,i)] = calcDifferenceandRatio(dbandbandwidth(:,i));
  [ddbandentropy(:,i), ddrbandentropy(:,i)] = calcDifferenceandRatio(dbandentropy(:,i));
  [ddbandenergy(:,i), ddrbandenergy(:,i)] = calcDifferenceandRatio(dbandenergy(:,i));
  [ddbandenergyratio(:,i), ddrbandenergyratio(:,i)] = calcDifferenceandRatio(dbandenergyratio(:,i));
  [ddbandflatness(:,i), ddrbandflatness(:,i)] = calcDifferenceandRatio(dbandflatness(:,i));
  [ddbandskewness(:,i), ddrbandskewness(:,i)] = calcDifferenceandRatio(dbandskewness(:,i));
  [ddbandflux(:,i), ddrbandflux(:,i)] = calcDifferenceandRatio(dbandflux(:,i));
  [ddbandspectraldecrease(:,i), ddrbandspectraldecrease(:,i)] = calcDifferenceandRatio(dbandspectraldecrease(:,i));
  
end


% create the band threshold vector (lower threshold from smoothed mean dB of spectrum)
%bandthresholds=ma+(md*bandsdthresh); % thresholded on global smoothed mean dB, not on dB of band!

% display thresholded Band (mean) Amplitudes
h=displayThresholdSpectogram(bandmean,freqbands,step,winLength,bandthresholds);
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click

if abs(saveimg)>0 % abs() as sometimes use -1 to only save method variant results, not spectral graphs
  export_fig('subband_ampl_thresh.jpg','-jpg',h)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sub band graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
% display Mean Band Amplitudes
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;
  plot(tsw(:),bandmean(:,i),'b-');
  %plot(tsw(:),bandstd(:,i),'r:');
  plot(tsw(:),bandmean(:,i)+bandstd(:,i),'r:');
  plot(tsw(:),bandmean(:,i)-bandstd(:,i),'r:');
  
  title(sprintf('Mean Amplitude & SD band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click



% display autocorrelation coeeficients for each band 
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  ylim([-1 1]);
  hold on;
  plot(tsw(:),bandacoeffs(:,i,1)/(max(bandacoeffs(:,i,1))-min(bandacoeffs(:,i,1))),'r-');
  % scale plots by 0th coefficient
  plot(tsw(:),bandacoeffs(:,i,2)/(max(bandacoeffs(:,i,1))-min(bandacoeffs(:,i,1))),'b-');
  plot(tsw(:),bandacoeffs(:,i,3)/(max(bandacoeffs(:,i,1))-min(bandacoeffs(:,i,1))),'g-');
  % scale each coeff by itself
  %plot(tsw(:),bandacoeffs(:,i,2)/(max(bandacoeffs(:,i,2))-min(bandacoeffs(:,i,2))),'b-');
  %plot(tsw(:),bandacoeffs(:,i,3)/(max(bandacoeffs(:,i,3))-min(bandacoeffs(:,i,3))),'g-');
  
  title(sprintf('Autocorr coeffs band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click

end

%%%%%%%%%%%%%%%%%%% display Band Features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if 0
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  ylim([0 1]);
  hold on;
  
  plot(tsw(:),bandenergy(:,i)/max(abs(bandenergy(:,i))),'k-');
  plot(tsw(:),bandbandwidth(:,i)/max(abs(bandbandwidth(:,i))),'g-');
  plot(tsw(:),bandentropy(:,i)/max(abs(bandentropy(:,i))),'b-');
  plot(tsw(:),bandflatness(:,i)/max(abs(bandflatness(:,i))),'m:');

  plot(tsw(:),bandacoeffs(:,i,1)/(max(bandacoeffs(:,i,1))-min(bandacoeffs(:,i,1))),'r-'); % autocorrel coef 0
  
  %plot(tsw(:),bandflux(:,i)/max(bandflux(:,i))-min(bandflux(:,i))),'c-');
  bfl=bandflux(:,i)/(max(bandflux(:,i))-min(bandflux(:,i))); % normalise
  bfl=bfl-min(bfl); % offset to positive
  plot(tsw(:),bfl,'c-');  
  
  title(sprintf('E, BW, entropy, flatness, autocorr 0, flux %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
h=gca;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click

if saveimg>0
  export_fig('subband_energy+BW+entropy+flat+autocorr0+flux.jpg','-jpg',h) % TODO: fix - this only saves 1 subband!!! save manually!
end
%end



%if 0
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  ylim([0 1]);
  hold on;

  plot(tsw(:),bandmean(:,i)/max(bandmean(:,i)),'b-');
  plot(tsw(:),bandenergyratio(:,i)/max(abs(bandenergyratio(:,i))),'r-');

  %plot(tsw(:),bandskewness(:,i)/(max(bandskewness(:,i))-min(bandskewness(:,i))),'g-');
  bsk=bandskewness(:,i)/(max(bandskewness(:,i))-min(bandskewness(:,i))); % normalise
  bsk=bsk-min(bsk); % offset to positive
  plot(tsw(:),bsk,'g-');
  
  %plot(tsw(:),bandspectraldecrease(:,i)/(max(bandspectraldecrease(:,i))-min(bandspectraldecrease(:,i))),'m-');
  bsd=bandspectraldecrease(:,i)/(max(bandspectraldecrease(:,i))-min(bandspectraldecrease(:,i))); % normalise
  bsd=bsd-min(bsd); % offset to 0-1
  plot(tsw(:),bsd,'m-');  

  title(sprintf('mean amplitude, ER, skewness, spdec %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
h=gca;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click

if saveimg>0
  export_fig('subband_mean+ER+skew+spdec.jpg','-jpg',h) % TODO: fix - this only saves 1 subband!!! save manually!
end
%end


%%%%%%%%%%%%%%%%%%%%%%%%% display Band Features in terms of SD from moving average %%%%%%%%%%%%%%%%%%%%%%

if 0
    
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  %ylim([0 1]);
  hold on;
  
  plot(tsw(:),scbandacoeffs(:,i),'b-');
  plot(tsw(:),thbandacoeffs(:,i)*max(scbandacoeffs(:,i)),'b*');

  title(sprintf('SD acoeff band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  %ylim([0 1]);
  hold on;
  
  plot(tsw(:),scbandenergy(:,i),'r-');
  plot(tsw(:),thbandenergy(:,i)*max(scbandenergy(:,i)),'r*');

  title(sprintf('SD Energy band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  %ylim([0 1]);
  hold on;
  
  plot(tsw(:),scbandenergyratio(:,i),'m-');
  plot(tsw(:),thbandenergyratio(:,i)*max(scbandenergyratio(:,i)),'m*');

  title(sprintf('SD ER band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  %ylim([0 1]);
  hold on;
  
  plot(tsw(:),scbandbandwidth(:,i),'g-');
  plot(tsw(:),thbandbandwidth(:,i)*max(scbandbandwidth(:,i)),'g*');

  title(sprintf('SD BW band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  %ylim([0 1]);
  hold on;
  
  plot(tsw(:),scbandentropy(:,i),'g-');
  plot(tsw(:),thbandentropy(:,i)*max(scbandentropy(:,i)),'g*');

  title(sprintf('SD Entropy band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  %ylim([0 1]);
  hold on;
  
  plot(tsw(:),scbandflatness(:,i),'b-');
  plot(tsw(:),thbandflatness(:,i)*max(scbandflatness(:,i)),'b*');

  title(sprintf('SD flatness band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  %ylim([0 1]);
  hold on;
  
  plot(tsw(:),scbandflux(:,i),'r-');
  plot(tsw(:),thbandflux(:,i)*max(scbandflux(:,i)),'r*');

  title(sprintf('SD flux band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  %ylim([0 1]);
  hold on;
  
  plot(tsw(:),scbandskewness(:,i),'m-');  
  plot(tsw(:),thbandskewness(:,i)*max(scbandskewness(:,i)),'m*');

  title(sprintf('SD skew band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  %ylim([0 1]);
  hold on;
  
  plot(tsw(:),scbandspectraldecrease(:,i),'g-');
  plot(tsw(:),thbandspectraldecrease(:,i)*max(scbandspectraldecrease(:,i)),'g*');

  %plot(tsw(:),scbandcentroid(:,i),'c-'); % not useful
  %plot(tsw(:),thbandcentroid(:,i)*max(scbandcentroid(:,i)),'c*');

  title(sprintf('SD specdec band %d - %d Hz',bandlo,bandhi));
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% graphs of vel & acc of each feature on all bands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0

% display velocity and acceleration of Mean Band Amplitudes
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;
  
  plot(tsw(:),[0; dbandmean(:,i)]/max(abs(dbandmean(:,i))),'rx:');
  plot(tsw(:),[0; 0; ddbandmean(:,i)]/max(abs(ddbandmean(:,i))),'c:');
  title(sprintf('Scaled Vel & Accel Mean Amplitude band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; dbandmean(:,i)],'rx:');
  %plot(tsw(:),[0; 0; ddbandmean(:,i)],'c:');
  %title(sprintf('Vel & Accel Mean Amplitude band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; drbandmean(:,i)]/max(abs(drbandmean(:,i))),'r:');
  %plot(tsw(:),[0; 0; ddrbandmean(:,i)]/max(abs(ddrbandmean(:,i))),'b:');
  %title(sprintf('Scaled Vel & Accel Ratio Mean Amplitude band %d - %d Hz',bandlo,bandhi));

  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click



% display velocity and acceleration of Bands' autocorrelation coefficient 0
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;
  
  plot(tsw(:),[0; dbandacoeffs(:,i)]/max(abs(dbandacoeffs(:,i))),'rx:');
  plot(tsw(:),[0; 0; ddbandacoeffs(:,i)]/max(abs(ddbandacoeffs(:,i))),'c:');
  title(sprintf('Scaled Vel & Accel Autocorrel band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; dbandacoeffs(:,i)],'rx:');
  %plot(tsw(:),[0; 0; ddbandacoeffs(:,i)],'c:');
  %title(sprintf('Vel & Accel Autocorrel band %d - %d Hz',bandlo,bandhi));

  %plot(tsw(:),[0; drbandacoeffs(:,i)]/max(abs(drbandacoeffs(:,i))),'rx:');
  %plot(tsw(:),[0; 0; ddrbandacoeffs(:,i)]/max(abs(ddrbandacoeffs(:,i))),'c:');
  %title(sprintf('Scaled Vel & Accel Ratio Autocorrel band %d - %d Hz',bandlo,bandhi));

  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click




% display velocity and acceleration of band centroid
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;
  
  plot(tsw(:),[0; dbandcentroid(:,i)]/max(abs(dbandcentroid(:,i))),'rx:');
  plot(tsw(:),[0; 0; ddbandcentroid(:,i)]/max(abs(ddbandcentroid(:,i))),'c:');
  title(sprintf('Scaled Vel & Accel Centroid band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; dbandcentroid(:,i)],'rx:');
  %plot(tsw(:),[0; 0; ddbandcentroid(:,i)],'c:');
  %title(sprintf('Vel & Accel Centroid band %d - %d Hz',bandlo,bandhi));

  %plot(tsw(:),[0; drbandcentroid(:,i)]/max(abs(drbandcentroid(:,i))),'r:');
  %plot(tsw(:),[0; 0; ddrbandcentroid(:,i)]/max(abs(ddrbandcentroid(:,i))),'b:');
  %title(sprintf('Scaled Vel & Accel Ratio Centroid band %d - %d Hz',bandlo,bandhi));

  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


% display velocity and acceleration of band bandwidth
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;
  
  plot(tsw(:),[0; dbandbandwidth(:,i)]/max(abs(dbandbandwidth(:,i))),'rx:');
  plot(tsw(:),[0; 0; ddbandbandwidth(:,i)]/max(abs(ddbandbandwidth(:,i))),'c:');
  title(sprintf('Scaled Vel & Accel Bandwidth band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; dbandbandwidth(:,i)],'rx:');
  %plot(tsw(:),[0; 0; ddbandbandwidth(:,i)],'c:');
  %title(sprintf('Vel & Accel Bandwidth band %d - %d Hz',bandlo,bandhi));

  %plot(tsw(:),[0; drbandbandwidth(:,i)]/max(abs(drbandbandwidth(:,i))),'r:');
  %plot(tsw(:),[0; 0; ddrbandbandwidth(:,i)]/max(abs(ddrbandbandwidth(:,i))),'b:');
  %title(sprintf('Scaled Vel & Accel Ratio Bandwidth band %d - %d Hz',bandlo,bandhi));

  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click



% display velocity and acceleration of band energy
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;
  
  plot(tsw(:),[0; dbandenergy(:,i)]/max(abs(dbandenergy(:,i))),'rx:');
  plot(tsw(:),[0; 0; ddbandenergy(:,i)]/max(abs(ddbandenergy(:,i))),'c:');
  title(sprintf('Scaled Vel & Accel Energy band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; dbandenergy(:,i)],'rx:');
  %plot(tsw(:),[0; 0; ddbandenergy(:,i)],'c:');
  %title(sprintf('Vel & Accel Energy band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; drbandenergy(:,i)]/max(abs(drbandenergy(:,i))),'r:');
  %plot(tsw(:),[0; 0; ddrbandenergy(:,i)]/max(abs(ddrbandenergy(:,i))),'b:');
  %title(sprintf('Scaled Vel & Accel Ratio Energy band %d - %d Hz',bandlo,bandhi));

  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click



% display velocity and acceleration of band energy ratio
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;
  
  plot(tsw(:),[0; dbandenergyratio(:,i)]/max(abs(dbandenergyratio(:,i))),'rx:');
  plot(tsw(:),[0; 0; ddbandenergyratio(:,i)]/max(abs(ddbandenergyratio(:,i))),'c:');
  title(sprintf('Scaled Vel & Accel energy ratio band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; dbandenergyratio(:,i)],'rx:');
  %plot(tsw(:),[0; 0; ddbandenergyratio(:,i)],'c:');  
  %title(sprintf('Vel & Accel energy ratio band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; drbandenergyratio(:,i)]/max(abs(drbandenergyratio(:,i))),'r:');
  %plot(tsw(:),[0; 0; ddrbandenergyratio(:,i)]/max(abs(ddrbandenergyratio(:,i))),'b:');
  %title(sprintf('Scaled Vel & Accel Ratio energy ratio band %d - %d Hz',bandlo,bandhi));

  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click



% display velocity and acceleration of band flatness
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;
  
  plot(tsw(:),[0; dbandflatness(:,i)]/max(abs(dbandflatness(:,i))),'rx:');
  plot(tsw(:),[0; 0; ddbandflatness(:,i)]/max(abs(ddbandflatness(:,i))),'c:');
  title(sprintf('Scaled Vel & Accel flatness band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; dbandflatness(:,i)],'rx:');
  %plot(tsw(:),[0; 0; ddbandflatness(:,i)],'c:');
  %title(sprintf('Vel & Accel flatness band %d - %d Hz',bandlo,bandhi));

  %plot(tsw(:),[0; drbandflatness(:,i)]/max(abs(drbandflatness(:,i))),'r:');
  %plot(tsw(:),[0; 0; ddrbandflatness(:,i)]/max(abs(ddrbandflatness(:,i))),'b:');
  %title(sprintf('Scaled Vel & Accel Ratio flatness band %d - %d Hz',bandlo,bandhi));
  
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click


% display velocity and acceleration of band skewness
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;

  title(sprintf('Scaled Vel & Accel skewness band %d - %d Hz',bandlo,bandhi));
  plot(tsw(:),[0; dbandskewness(:,i)]/max(abs(dbandskewness(:,i))),'rx:');
  plot(tsw(:),[0; 0; ddbandskewness(:,i)]/max(abs(ddbandskewness(:,i))),'c:');
  
  %plot(tsw(:),[0; dbandskewness(:,i)],'rx:');
  %plot(tsw(:),[0; 0; ddbandskewness(:,i)],'c:');
  %title(sprintf('Vel & Accel skewness band %d - %d Hz',bandlo,bandhi));

  %plot(tsw(:),[0; drbandskewness(:,i)]/max(abs(drbandskewness(:,i))),'r:');
  %plot(tsw(:),[0; 0; ddrbandskewness(:,i)]/max(abs(ddrbandskewness(:,i))),'b:');
  %title(sprintf('Scaled Vel & Accel Ratio skewness band %d - %d Hz',bandlo,bandhi));

  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click



% display velocity and acceleration of band flux
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;
  
  plot(tsw(:),[0; dbandflux(:,i)]/max(abs(dbandflux(:,i))),'rx:');
  plot(tsw(:),[0; 0; ddbandflux(:,i)]/max(abs(ddbandflux(:,i))),'c:');
  title(sprintf('Scaled Vel & Accel flux band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; dbandflux(:,i)],'rx:');
  %plot(tsw(:),[0; 0; ddbandflux(:,i)],'c:');
  %title(sprintf('Vel & Accel flux band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; drbandflux(:,i)]/max(abs(drbandflux(:,i))),'r:');
  %plot(tsw(:),[0; 0; ddrbandflux(:,i)]/max(abs(ddrbandflux(:,i))),'b:');
  %title(sprintf('Scaled Vel & Accel Ratio flux band %d - %d Hz',bandlo,bandhi));
  
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click



% display velocity and acceleration of band spectral decrease
clf reset % clear graph
bandhi=-1;
for i=1:length(freqbands)
  bandlo=bandhi+1;
  bandhi=freqbands(i);
  subplot(floor((length(freqbands)+1)/2), 2, i); % tile plots with each plot taking up a row
  hold on;
  
  plot(tsw(:),[0; dbandspectraldecrease(:,i)]/max(abs(dbandspectraldecrease(:,i))),'rx:');
  plot(tsw(:),[0; 0; ddbandspectraldecrease(:,i)]/max(abs(ddbandspectraldecrease(:,i))),'c:');
  title(sprintf('Scaled Vel & Accel spectral decrease band %d - %d Hz',bandlo,bandhi));
  
  %plot(tsw(:),[0; dbandspectraldecrease(:,i)],'rx:');
  %plot(tsw(:),[0; 0; ddbandspectraldecrease(:,i)],'c:');
  %title(sprintf('Vel & Accel spectral decrease band %d - %d Hz',bandlo,bandhi));

  %plot(tsw(:),[0; drbandspectraldecrease(:,i)]/max(abs(drbandspectraldecrease(:,i))),'r:');
  %plot(tsw(:),[0; 0; ddrbandspectraldecrease(:,i)]/max(abs(ddrbandspectraldecrease(:,i))),'b:');
  %title(sprintf('Scaled Vel & Accel Ratio spectral decrease band %d - %d Hz',bandlo,bandhi));
  
  %xlabel('time (sec)');
  %ylabel('Amplitude');
  grid on;
end
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click

end


%%%%%%%%%%%%%%%%%% peform heuristics for determining if the sound is an approaching vehicle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for approaching vehicle energy, bandwidth, entropy, auto correlation should be high
%                         flatness, flux, abs(spectral decrease) should be low
%  many of these features are interrelated - from smoothed (1 sec window) graphs it appears that
%  ignoring offsets:
%   - flux follows flatness but is less detailed (flatter with just a few peaks/valleys)
%   - spectral decrease mirrors flatness
%   - bandwidth follows spectral decrease (but is flatter)
%   - entropy follows bandwidth
%   - entropy also follows energy to a large degree
%
%  ie energy -> entropy -> bandwidth -> flatness -> mirrored spectral decrease & flattened flux
%
% so, as signal amplitude (energy) is already utilised in a more focused manner (sub bands), these
% global features add no additional information
% however, it may be that these features just need to be viewed in another way to make them independently useful

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% count the number (N) of "active" bands (those with mean of freq components above the threshold) in bands 2-11 (40-1000Hz)
% and if band 11 is active, then also include bands 12-14 (1-2,2-3 kHz).
% also note the top active band (F) upto 14
% if N >= 4 (for example) ie if 40% of the bands are active 
%           AND (F >= 9 (ie top active freq band > 301-400Hz)
%               OR band 1 is inactive (then prob not wind))
% then pass to next stage to determine if it's a vehicle
% note - that in a real time-based implementation, will want to see if N is
% increasing (if low value) or just a spike (even if high value) to warn of approach


activebands=zeros(size(yd,1), length(freqbands));
highestactiveband=zeros(size(yd,1),1);

for i=1:length(freqbands)
   activebands(:,i)=bandmean(:,i)>=bandthresholds(:,i);
   a=find(activebands(:,i)>0); % indexes (in time) of when band is active
   if i<15 % only count bands < 3kHZ
     highestactiveband(a)=i;
   end
end

vact=sum(activebands(:,2:11),2); % sum of active bands 2-11
vact=vact+(activebands(:,11).*sum(activebands(:,12:14),2)); % if band 11 is active, then adds bands 12-14
vact=vact.*(highestactiveband>=9 | activebands(:,1)==0); % active if top active >=9 OR band 1 inactive

vact=min(vact,10); % only count up to 10
%vact(find(vactl<4))=0; % ignore activity counts less than 4

% this signal can have high amplitude spikes which may cause false reporting of vehicle approach
% so these spikes should be eliminated

% eliminate (noise) spikes using a median filter
%vactbampl = medfilt1(vact,medfiltw);

% eliminate spikes in activity caused by noise ie spikes eminating from 0 and then dropping back to zero
vactbampl=removeSpikes(vact,minvactw,0); % affects all spikes


if 0

clf reset % clear graph
hold on
title('Vehicle Likelyhood based on Band Amplitude above Normal');
%title('Normalised Signal Power, Velocity and Acceleration');
xlabel('time (sec)');
ylabel('Count');
ylim([0 100]);

plot(tsw(:),vactbampl(:)*10,'r-'); % scale to 0-100
plot(tsw(:),acoeffs(:,1)/(max(acoeffs(:,1))-min(acoeffs(:,1)))*100,'b-.'); % auto correlation coefficient 0 (scale 0-100)
%plot(tsw(:),rms(:)/max(rms)*100,'g:'); % this represents signal amplitude (scale 0-100)
lgn=legend('vehicle activity','autocorr coeff',4);%'RMS',4);
grid on;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test for approaching vehicle using thresholded band spectral features (energy, bandwidth, flatness + autocorr 0 & flux for imminence)
% based on SD away from moving average
%
% for approaching vehicle - energy and bandwidth should be above SD for most bands
% in 401-2000Hz (bands 9-13) with bands 250-400Hz (bands 7,8), 2kHz - 11kHz (bands 14-16) as optional:
%    to consider as vehicle, either 4 out of 5 main bands, or 3 out of 5 main bands plus 3 out of 5 optional bands
% and flatness (scaled to mean and SD) should be < 0 for most (4 out of 6) bands >= 750Hz (bands 11-16)

vae=zeros(size(yd,1), 1);
vaeopt=zeros(size(yd,1), 1);
vab=zeros(size(yd,1), 1);
vabopt=zeros(size(yd,1), 1);
vaf=zeros(size(yd,1), 1);
ivaa=zeros(size(yd,1), 1);
ivaf=zeros(size(yd,1), 1);

for i=9:13 % bands covering 401Hz - 2Khz
  vae=vae+(scbandenergy(:,i)>=bandsdthresh*rofactor); % increment count for energy significantly high in band
  vab=vab+(scbandbandwidth(:,i)>=bandsdthresh*rofactor); % increment count for bandwidth significantly high in band
end

for i=[7 8 14 15 16] % bands covering 250-400Hz and 2Khz-11kHz
  vaeopt=vaeopt+(scbandenergy(:,i)>=bandsdthresh*rofactor); % increment count for energy significantly high in band
  vabopt=vabopt+(scbandbandwidth(:,i)>=bandsdthresh*rofactor); % increment count for bandwidth significantly high in band
end

for i=11:16 % bands with freq > 750Hz
  vaf=vaf+(scbandflatness(:,i)<0); % increment count for flatness less than mean (ie reasonably low)
end

vacte=vae>=4 | (vae==3 & vaeopt>=3); % 4 out of 5 main bands, or 3 out of 5 main bands and most optional bands
vactb=vab>=4 | (vab==3 & vabopt>=3); % 4 out of 5 main bands, or 3 out of 5 main bands and most optional bands
vactf=vaf>=4; % more than half (actually at least 2/3rds) the important bands

%vact = vacte & vactb & vactf; % all the above conditions must be true
vact = vacte + vactb + vactf; % add the 3 features
vact=removeSpikes(vact,minvactw,0); % remove all spikes!

% for IMMINENT vehicle, add autocorr above SD for most bands >= 251Hz, and
% flux low for most bands 251Hz - 2kHz

for i=[7:16] % bands covering 251-11kHz
  ivaa=ivaa+(scbandacoeffs(:,i)>=bandsdthresh*rofactor); % increment count for auto correl significantly high in band
end

for i=[7:12] % bands covering 251-2kHz
  ivaf=ivaf+(scbandflux(:,i)<=-bandsdthresh*rofactor); % increment count for flux significantly low in band
end

ivacta=ivaa>5; % at least 6 out of 10 bands above normal
ivactf=ivaf>3; % at least 4 out of 6 bands below normal

% values >= 3 should indicate a vehicle approaching, with values >3 indicating immininence!
% not worrying about additional filtering out of spikes for imminence as these are only added if main 3 features are active
vactbfeat=vact+((vact>0).*(ivacta+ivactf)); % add the immanence detectors if all the main features are positive

vactbfeat=vactbfeat*2; % scale from 0-5 to 0-10 to make equal importance to vactbampl
% eliminate spikes in activity caused by noise ie spikes eminating from 0 and then dropping back to zero
% this doesn't really appear to be necessary here, it is just a precaution
vactbfeat=removeSpikes(vactbfeat,minvactw,0);

% both band amplitude and spectral features should be present at the same
% however, it seems that thresholding based on roll off affects band amplitude more than
% features, so, if amplitude is 0, halve the feature score - leaving it
% still intact - as it it the first thing to register approaching vehicle
% but if only amplitude is active - leave it as it has fewer false positives

vscore=(vactbfeat(:)+((vactbfeat(:)).*(vactbampl>0)))/2;
vscore=vscore+vactbampl(:);
vscore=vscore/2; % scale score to 0-10

%((vactbampl(:)+vactbfeat(:))+((vactbampl(:)+vactbfeat(:)).*(vactbampl>0 & vactbfeat>0)))/4;
%vscore=((vactbampl(:)+vactbfeat(:)).*(vactbampl>0 & vactbfeat>0))/2; % scale score to 0-10
%vscore=(vactbampl+vactbfeat)/2; % scale score to 0-10

clf reset % clear graph
hold on
title('Vehicle Likelyhood based on Band Amplitude & Features over Thresholds');
xlabel('time (sec)');
ylabel('Count');
ylim([0 100]);

plot(tsw(:),50+(vactbampl(:)*5),'b-'); % scale to 0-50 and put at top of graph
plot(tsw(:),0+(vactbfeat(:)*5),'g-'); % scale to 0-50 and put at bottom of graph
plot(tsw(:),vscore(:)*10,'r-'); % scale score to 0-100

plot(tsw(:),rofactor(:)*10,'m.-'); % scale score to 0-100
lgn=legend('band ampl','band feat','score','rofactor x 10',4);%'RMS',4);


%plot(tsw(:),acoeffs(:,1)/(max(acoeffs(:,1))-min(acoeffs(:,1)))*10,'m-'); % auto correlation coefficient 0 (scale 0-10)
%lgn=legend('band ampl','band feat','score','autocorr coeff',4);%'RMS',4);

%plot(tsw(:),rms(:)/max(rms)*100,'c:'); % this represents signal amplitude (scale 0-100)

grid on;
h=gca;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
if abs(saveimg)>0 % abs() as sometimes use -1 to only save results and not spectral graphs
  export_fig('feat+ampl_score.jpg','-jpg',h)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for recognition testing experiments, rather than actual application,
% scores are done slightly differently:
% - both band amplitude and spectral features should be present at the same
%   if they aren't, then halve the score

vehfscore=(vactbfeat(:)+((vactbfeat(:)).*(vactbampl>0)))/2; % 0-10
vehascore=(vactbampl(:)+((vactbampl(:)).*(vactbfeat>0)))/2; % 0-10
vehscore=(vehfscore+vehascore)*5; % scale score to 0-100

% remove all spikes and drop outs!
vehscore=removeSpikes(vehscore,minvactresw,0);
vehscore=removeDropOuts(vehscore,minvactresw);


% score 0 - no detection: below 10% (either nothing or some partial feat only or low mag only)
% score 1 - low level detection: 10-29% (feat only is 0-25, mag only is 0-25)
% score 2 - medium level detection: 30-59% (both feat and mag present, but low mag and no immanence)
% score 3 - high level detection: 60%-100% (both present with good strength, either immanence triggered
%                                 and or magnitude high)

% raise each detection event to its maximum
peakedvehscore=raiseToPeak(vehscore);

% score each detection event 0-3 based on its maximum
detectionscore=(peakedvehscore(:)>=10)+(peakedvehscore(:)>=30)+(peakedvehscore(:)>=60);


% graph detection score

clf reset % clear graph
hold on
title('Vehicle Detection Score (low, med, high)');
xlabel('time (sec)');
ylabel('score');
ylim([-1 3]);

plot(tsw(:),detectionscore(:),'r-');
plot(ts(:),x(:)/max(abs(x)),'g:');

lgn=legend('detection score','signal',4);

grid on;
h=gca;
disp('click graph to continue')
k = waitforbuttonpress(); % pause on keystroke or mouse click
if abs(saveimg)>0 % abs() as sometimes use -1 to only save results and not spectral graphs
  export_fig('detection_score.jpg','-jpg',h)
end

% count the total number of high, medium, and low vehicle detection events
[vhicnt,vhiidx,vhiwidth]=findPeaks(peakedvehscore(:)>=60,0);
[vmedcnt,vmedidx,vmedwidth]=findPeaks((peakedvehscore(:)<60)&(peakedvehscore(:)>=30),0);
[vlocnt,vloidx,vlowidth]=findPeaks((peakedvehscore(:)<30)&(peakedvehscore(:)>=10),0);


sprintf('results for sound file: %s',fileName)

for i=1:vhicnt
  sprintf('significant vehicle %d detected at %.2f secs, lasting %.2f secs',i,vhiidx(i)*step,vhiwidth(i)*step)
end

for i=1:vmedcnt
  sprintf('other vehicle %d detected at %.2f secs, lasting %.2f secs',i,vmedidx(i)*step,vmedwidth(i)*step)
end

for i=1:vlocnt
  sprintf('possible vehicle %d detected at %.2f secs, lasting %.2f secs',i,vloidx(i)*step,vlowidth(i)*step)
end

