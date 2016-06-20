   %% Time specifications:
   
   %N og x
   
   Fs = 1/TR;                      % samples per second
   dt = 1/Fs;                     % seconds per sample
   t = (0:dt:N*dt-dt)';
   

   %% Fourier Transform:
   Xf = fftshift(fft(x));
   
   %% Frequency specifications:
   dF = Fs/N;                      % hertz
   f = 0:dF:Fs/2-dF;           % hertz
   
   %% Plot the spectrum:
  
   plot(f,abs(Xf(end/2+1:end))/N);
%    xlabel('Frequency (in hertz)');
%    title('Magnitude Response');