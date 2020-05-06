function [ dl, dt, L, time_fieldsData, complx_freq_data, freq_bin] = Compute_FFT_on_timeDomain_Data( filename , range)

    start = 1;
    begin_pos = 3;

    time_fieldsData = importdata(filename,' '); 
    
    dl = time_fieldsData(1:1,1);         % mesh discretization length
    num_of_iter = time_fieldsData(2:2,1);   % total number of samples / num_of_iter
    if (range< num_of_iter)
        num_of_iter = range;
    end
    time_fieldsData = time_fieldsData(begin_pos:num_of_iter,1);
    
    num_of_iter = length(time_fieldsData);
    
    L = num_of_iter;
    dl = dl*1e-3;
    c = 299792458;                   % speed of light
    dt = 0.5*dl/c;                   % maximum time step for SCN
%     dt = dl/(sqrt(2)*c);             % maximum time step for 2D simulations

    % Computing the FFT 
    fs = 1/dt; 
    N = L;
    f_res = fs/N;
    f_max = f_res * (N-1); 
    freq_bin = 0:f_res:f_max; % frequency bins
    freq_bin = freq_bin.';

    winvec = 1; %hamming(N); %windowing function
    complx_freq_data = fft( time_fieldsData.*winvec );
    
end % function end
% 
% figure (1)
% plot(ffG,10*log10(Ey1),'g')
% hold on
% title('Amplitude versus Frequency  ');
% xlabel('Frequency GHz');
% ylabel('Amplitude');
% 
% % 
% % % figure(3) 
% % % plot(t(start:end)/dt,10*log10(abs(EY1(start:end))),'g',...
% % %      t(start:end)/dt,10*log10(abs(EY2(start:end))) ,'b',...
% % %      t(start:end)/dt, 10*log10(abs(EY3(start:end))),'r')
% % %  
% % %  hold on 
% % % title(', ');
% % % xlabel('Time steps');
% % % ylabel('Ey');
% % 
% % % dl=dl*1e3;
% % ref_max = max( abs(EY1) );
% % db_max = max( abs(EY2) );
% % mpml_max = max( (EY3) );
% % 
% % % ref_max = 1; %max( EY1 );
% % % db_max = 1; %max( EY2 );
% % % mpml_max = 1; %max( EY3 );
% % 
% % figure(fig2) 
% % plot(t(start:end)/dt,EY1(start:end)/ref_max,'g',...
% %      t(start:end)/dt,EY2(start:end)/db_max ,'b',...
% %      t(start:end)/dt, EY3(start:end)/mpml_max,'r')
% %  
% %  hold on 
% % title(', ');
% % xlabel('Time steps');
% % ylabel('Ey');
% % 
% % 
% % 
% % % 
% % % % MORE PLOTTING OPTIONS - MULTI PLOTs etc.
% % % 
% % figure(fig1);
% % subplot(4,1,1);
% % plot(   t(start:end)/dt,   EY1(start:end)/ref_max,   'g') 
% % title (' Normalized Field Values ');
% % xlabel('Time step ');
% % ylabel('V/m');
% % % 
% % % subplot(4,1,2); 
% % % plot(     t(start:end)/dt,     EY2(start:end)/db_max ,    'b') 
% % % xlabel('Time step'); ylabel('V/m');
% % % 
% % % subplot(4,1,3);
% % % plot(  t(start:end)/dt,    EY3(start:end)/mpml_max,   'r') 
% % % xlabel('Time step'); ylabel('V/m');
% % % 
% % % subplot(4,1,4);
% % % plot(  t(start:end)/dt,    EY1(start:end)/ref_max,   'g') 
% % % xlabel('Time step'); ylabel('V/m');
% % % hold on 
% % % plot(     t(start:end)/dt,     EY2(start:end)/db_max ,    'b') 
% % % plot(     t(start:end)/dt,     EY3(start:end)/mpml_max ,    'r') 
% % % 
% % % 
% % % % % 
% % % figure(fig2);
% % % subplot(3,1,1);
% % % plot(   t(start:end)/dt,   EY1(start:end),   'g') 
% % % xlabel('Time step ');
% % % ylabel('V/m');
% % % 
% % % subplot(3,1,2); 
% % % plot(     t(start:end)/dt,     EY2(start:end) ,    'b') 
% % % xlabel('Time step'); ylabel('/m');
% % % 
% % % subplot(3,1,3);
% % % plot(  t(start:end)/dt,    EY3(start:end),   'r') 
% % % xlabel('Time step'); ylabel('V/m');
% % % % 
% % 
% % [E,tp] = findpeaks(20*log10(abs(EY1-EY2) ));
% % Error = 20*log10(((abs(EY3-EY1)/mpml_max )));
% % % Error = movmean(Error,4);
% % 
% % % % % % 
% % figure(fig3) 
% % plot(t(1:end)/dt,Error)
% % hold on
% % title('difference');
% % xlabel('Time steps');
% % % ylim([-200,-50])
% % % xlim([0,220])
% % 
% % % 
% % % figure(12)
% % % plot(tp,E)
% % % hold on
% % % title('difference');
% % % xlabel('Time steps');
% % % % % ylabel('Is');
% 
% end
