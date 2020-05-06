% This code is part of the GEMINI package

% Author: J.Odeyemi
% Email: Jomiloju.odeyemi@nottingham.ac.uk

% ------------------------------------------------------------------------------------------------------------
% This script plots time  and frequency domain samples from multiple data files.
% 
%------------------------------------------------------------------------------------------------------------

c = 3e8;
bandwidth = 1e9;
centre_freq = 1e9;
dl = 1e-3;
dt = dl/(c*2);
num_of_iter = 1e4;

[ xt , xf, xf_bins] = Gaussian_Excitation( dt, num_of_iter, bandwidth, centre_freq );

ampl_fft = abs(xf)/num_of_iter;

figure(2);
subplot(2,1,1); 
plot(0:num_of_iter, xt,'b'); 
title('Time Domain Plot of Excitation');
xlabel('Time step'); 
ylabel('Amplitude');
   

subplot(2,1,2); 
plot(xf_bins, ampl_fft,'r'); 
title('Frequency Domain Plot of Excitation ');      
xlabel('Frequency (Hz)'); 
ylabel('Magnitude |X(f)|');