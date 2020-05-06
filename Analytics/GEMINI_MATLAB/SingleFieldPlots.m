% This code is part of the GEMINI package
% Author: J.Odeyemi
% Email: Jomiloju.odeyemi@nottingham.ac.uk
%
%------------------------------------------------------------------------------------------------------------
% This script plots the time and frequency domain from a single data file
%
%------------------------------------------------------------------------------------------------------------

filename  = 'electric_field_rectangular_WG.txt';
range = 20000;
[ dl, dt, num_of_iter, time_fieldsData, complx_freq_data, freq_bin] = Compute_FFT_on_timeDomain_Data(filename, range);

%------------------------------------------------------------------------------------------------------------

ampl_fft = abs(complx_freq_data)/num_of_iter; % amplitude of complex freq data

figure(2);
subplot(2,1,1); 
plot(1:num_of_iter, time_fieldsData,'b'); 
title('Time Domain Field Plot');
xlabel('Time step'); 
ylabel('Amplitude');
   

subplot(2,1,2); 
plot(freq_bin(1:num_of_iter/2), ampl_fft(1:num_of_iter/2),'r'); 
title('Magnitude of FFT');      
xlabel('Frequency (Hz)'); 
ylabel('Magnitude |X(f)|');

%------------------------------------------------------------------------------------------------------------

