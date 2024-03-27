function CreateStudyStruc_BA(mattype,filename,outputdir,data,na,np,Nx,Ny,shift,lw)
    if na==100
        load('studytofill_denoising.mat');
    end 
    
    study.path=outputdir; %destination folder
    study.filename=filename;
    study.liststring=[study.path '\' study.filename];
    study.params.np = np;
    study.multiplicity = na;
    study.shift = shift;
    study.lw = lw;

    study.data.real=zeros(study.multiplicity,Nx,Ny,study.params.np); %real part of the data in line
    study.data.imag=zeros(study.multiplicity,Nx,Ny,study.params.np); %real part of the data in line

    %fill in your data here
    if mattype=="noiseless"
        data2=repmat(data,na,1);
        study.data.real =real(data2);
        study.data.imag =imag(data2);
    elseif mattype=="noisy"
        study.data.real =real(data);
        study.data.imag =imag(data);
    end 
        


    t=datetime('now');
    study.time=datestr(t);

    if(~exist(outputdir,"dir"))
        mkdir(outputdir);
    end
    save([outputdir '/' study.filename '.mat'],'study','-v7.3')


end

