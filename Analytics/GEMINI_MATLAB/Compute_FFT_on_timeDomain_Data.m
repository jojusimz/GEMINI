
% This code is part of the GEMINI package
% Author: J.Odeyemi
% Email: Jomiloju.odeyemi@nottingham.ac.uk
%
%------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------

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
