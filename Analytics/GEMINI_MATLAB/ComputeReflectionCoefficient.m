

%total_field = 'electric_field_parallel_plateWG.txt';
%reference_field = 'reference_electric_field_parallel_plateWG.txt';

total_field = 'electric_field_parallel_plateWG2.txt';
reference_field = 'reference_electric_field_parallel_plateWG2.txt';
range  = 900;

%--------------------------------------------------------------------------------------------------------------------

[num_of_iter, reflec_xt, reflec_xf, freq_bin] = Compute_Reflection_Coefficient(total_field, reference_field, range);

ampl_of_reflect = abs(reflec_xf)/num_of_iter;

%--------------------------------------------------------------------------------------------------------------------

figure(2);
subplot(2,1,1); 
plot(1:num_of_iter, reflec_xt,'b'); 
title('Time Domain Numerical reflection coefficient extracted from a TLM simulation');
xlabel('Time step'); 
ylabel('Amplitude');
hold on

subplot(2,1,2); 
plot(freq_bin(1:num_of_iter/2),20*log10(ampl_of_reflect(1:num_of_iter/2)));
title('Frequency Domain Numerical reflection coefficient extracted from a TLM simulation ');
xlabel('Frequency (GHz)');
ylabel('Reflection Coefficient (dB)');
xlim([20e9,40e9])
hold on
   
