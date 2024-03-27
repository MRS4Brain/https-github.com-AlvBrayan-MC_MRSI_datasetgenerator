function  Rephased_Data_t  = Compute_MinPhase( Data_t)
%COMPUTE_MINPHASE Summary of this function goes here

abs_Data_f=abs(fft(Data_t));
Rephased_Data_t=ifft(abs_Data_f .* exp (-1j * imag (hilbert (log(abs_Data_f))))); 


end

