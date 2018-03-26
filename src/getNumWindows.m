function nWindows = getNumWindows(x,winStep,winLength)

  nWindows=floor((length(x)-winLength)/winStep) + 1;

