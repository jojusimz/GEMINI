% This code is part of the GEMINI package
% Author: J.Odeyemi
% Email: Jomiloju.odeyemi@nottingham.ac.uk
%
%------------------------------------------------------------------------------------------------------------
% This script plots the time and frequency domain from a single data file
%
%------------------------------------------------------------------------------------------------------------

function [num_of_iter1, reflec_xt, reflec_xf, ref_freq_bin] = Compute_Reflection_Coefficient(filename1, filename2, range)

%total field data
    [ dl1, dt1, num_of_iter1, ttl_time_fieldsData, ttl_complx_freq_data, ttl_freq_bin] = Compute_FFT_on_timeDomain_Data(filename1, range);

%reference field data
    [ dl2, dt2, num_of_iter2, ref_time_fieldsData, ref_complx_freq_data, ref_freq_bin] = Compute_FFT_on_timeDomain_Data(filename2, range );

% reflected in time 
    reflec_xt = ttl_time_fieldsData - ref_time_fieldsData;               

% reflected in frequency
    reflec_xf =(ttl_complx_freq_data- ref_complx_freq_data)./ref_complx_freq_data;      % 
    
end

%     % %PLOTS
%     figure (1)
%     plot(ff/1e9,20*log10(R));
%     hold on
%     % title('Numerical reflection coefficient extracted from a TLM simulation of FW centred-fed dipole antenna situated in computational domain of varying sizes with matched termination ');
%     title('Reflection coefficient versus Frequency');
%     xlabel('Frequency (GHz)');
%     ylabel('Reflection Coefficient (dB)');
%     xlim([20 42])
%     % ylim([-650 -200]); 
% 
%     figure (3)
%     plot(t/dt,EY_total-EY_ref);
%     hold on
%     title('Reflected field');
%     xlabel('time');
%     ylabel('EY');
%     % 
%     [E,tp] = findpeaks(20*log10(abs(EY_total-EY_ref) / max(abs(EY_ref))));
%     Error = 20*log10(((abs(EY_total-EY_ref) / max(abs(EY_ref)))));
% 
% 
% % [E,tp] = findpeaks (E);
% 
% % figure (223)
% % plot(  ff/1e9, R_arg/pi);
% % hold on
% % title('Phase');
% % xlabel('Frequency');
% % ylabel('Phase');
% % xlim([20,40])
% % 
% % figure (225)
% % plot(  t/dt, Error);
% % hold on
% % title('Reflected wave Ey');
% % xlabel('time');
% % ylabel('EY');
% % % xlim([0 1000])
% % % ylim([-200 -20])
% % 
% % figure (226)
% % plot(t(1:end)/dt,EY_total(1:end));
% % hold on
% % title('Total field');
% % xlabel('time');
% % ylabel('EY');
% 
% % 
% % 
% % 
% % %FFT 2
% % % E_back_fft = fft( E_back.*winvec );
% % % E_back_fft = (abs( E_back_fft )/ L);
% % % 
% % % E_forward_fft = fft( Ey_ref.*winvec );
% % % E_forward_fft = (abs( E_forward_fft )/ L);
% % 
% % % %PLOTS
% % % figure (5)
% % % plot(ff,10*log10(E_back_fft));
% % % % hold on
% % % title('FFT of backward wave for dl = 0.6mm ');
% % % xlabel('Frequency Hz');
% % % ylabel('Ey in dB');
% % 
% % % %PLOTS
% % % figure (6)
% % % plot(ff,10*log10(E_forward_fft));
% % % % hold on
% % % title('FFT of forward wave for dl = 0.6mm ');
% % % xlabel('Frequency Hz');
% % % ylabel('Ey in dB');
% % 
% % 
% % 
% % 
% %PLOTS
% figure (4)
% plot(t(1:end)/dt,r);
% % % hold on
% title('reflection coefficient for dl = 0.6mm ');
% xlabel('Frequency Hz');
% ylabel('reflection coefficient');
