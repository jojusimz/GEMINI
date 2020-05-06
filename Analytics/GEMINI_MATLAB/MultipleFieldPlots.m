% This code is part of the GEMINI package

% Author: J.Odeyemi
% Email: Jomiloju.odeyemi@nottingham.ac.uk

%------------------------------------------------------------------------------------------------------------
% This script plots time  and frequency domain samples from multiple data files.
% 
%------------------------------------------------------------------------------------------------------------

filename1  = 'electric_field1.txt';
filename2  = 'electric_field2.txt';
filename3  = 'electric_field3.txt';

range = 10000000;

[ dl1, dt1, num_of_iter1, time_fieldsData1, complx_freq_data1, freq_bin1] = Compute_FFT_on_timeDomain_Data(filename1, range);
[ dl2, dt2, num_of_iter2, time_fieldsData2, complx_freq_data2, freq_bin2] = Compute_FFT_on_timeDomain_Data(filename2, range);
[ dl3, dt3, num_of_iter3, time_fieldsData3, complx_freq_data3, freq_bin3] = Compute_FFT_on_timeDomain_Data(filename3, range);

%------------------------------------------------------------------------------------------------------------

ampl_fft1 = abs(complx_freq_data1)/num_of_iter1; % amplitude of complex freq data
ampl_fft2 = abs(complx_freq_data2)/num_of_iter2; % amplitude of complex freq data
ampl_fft3 = abs(complx_freq_data3)/num_of_iter3; % amplitude of complex freq data
%------------------------------------------------------------------------------------------------------------

figure(2);
subplot(2,1,1); 

plot(1:num_of_iter1, time_fieldsData1,'g',...
     1:num_of_iter2, time_fieldsData2,'b',...
     1:num_of_iter3, time_fieldsData3,'r')
title('Time Domain Field Plot');
xlabel('Time step'); 
ylabel('Amplitude');
   

subplot(2,1,2); 
plot(freq_bin1(1:num_of_iter1/2), ampl_fft1(1:num_of_iter1/2),'g',...
     freq_bin2(1:num_of_iter1/2), ampl_fft2(1:num_of_iter1/2),'b',...
     freq_bin3(1:num_of_iter1/2), ampl_fft3(1:num_of_iter1/2),'r')

title('Magnitude of FFT');      
xlabel('Frequency (Hz)'); 
ylabel('Magnitude |X(f)|');

%------------------------------------------------------------------------------------------------------------

