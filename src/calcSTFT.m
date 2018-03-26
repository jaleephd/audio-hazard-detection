function y = calcSTFT(x,winStep,winLength)

% create Hamming window
hWin = hamming(winLength);

nWindows=getNumWindows(x,winStep,winLength);

iter=0;

if iter>0 % iterate through each time slice - slow, just do as eg for clarity
  
  % perform an FFT on each window
  for i=1:nWindows
    xoffset = ((i-1)*winStep)+1;
    XX=x(xoffset:xoffset+winLength-1).*hWin; % apply a Hamming window to this time slice
    % perform DFT on window
    YY = fft(XX); % contains + and -ve freq
    % extract +ve freq components
    if (mod(winLength,2)==1)
      k=(winLength+1)/2;
    else
      k=winLength/2;
    end
    YY=YY(1:k,:); % keep the +ve components
    y(i,:)=YY; % and store for this window
  end

else % do FFT on entire matrix in one go (faster)

  % build matrix of windowed slices
  XX=zeros(winLength,nWindows);
  for i=1:nWindows
    xoffset = ((i-1)*winStep)+1;
    XX(:,i) = x(xoffset:xoffset+winLength-1).*hWin; % apply a Hamming window
  end
  % compute fast fourier transform on all the slices
  YY = fft(XX);

  % extract +ve freq components
  if (mod(winLength,2)==1)
    k=(winLength+1)/2;
  else
    k=winLength/2;
  end
  y=YY(1:k,:)';

end


