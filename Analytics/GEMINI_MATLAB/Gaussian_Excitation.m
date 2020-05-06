% This code is part of the GEMINI package
% Author: J.Odeyemi
% Email: Jomiloju.odeyemi@nottingham.ac.uk
%
%------------------------------------------------------------------------------------------------------------
% This script plots the time and frequency domain from a single data file
%
%------------------------------------------------------------------------------------------------------------



function [ xt , xf, xf_bins] = Gaussian_Excitation( dt, num_of_iter, bandwidth, centre_freq )

    %sampling rate is determined by 1/dt.
   
    fs=1/dt; %sampling frequency
    fn=fs/2;
    t=(0:num_of_iter); %time base
    xt = (dt*t-(1*dt)).*bandwidth.* exp(-(bandwidth*(dt*t-(10*dt))).^2).*sin(2*pi*centre_freq*dt*t);

    xf =fftshift(fft(xt));
    xf_bins = fs*(-num_of_iter/2:num_of_iter/2)/num_of_iter; %Frequency Vector
    
end
