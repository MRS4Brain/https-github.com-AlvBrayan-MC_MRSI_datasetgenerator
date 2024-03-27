function [ Reph_CData_t ] = Rephase(  CData_t  )
    CData_f=fft(CData_t);
    pha=smooth(unwrap(angle(CData_f)),100);
    Reph_CData_f=CData_f(:).*exp(-1i*pha);
    
    pha=smooth(unwrap(angle(Reph_CData_f)),100);
    Reph_CData_f=Reph_CData_f(:).*exp(-1i*pha);
    
    pha=smooth(unwrap(angle(Reph_CData_f)),100);
    Reph_CData_f=Reph_CData_f(:).*exp(-1i*pha);
    
    Reph_CData_t=ifft(Reph_CData_f);
end

