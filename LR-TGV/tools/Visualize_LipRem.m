function Visualize_LipRem(name_fig,Mask2dBrain,MinMax_pt,legends,Vols_tkk )

%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here


nDims=ndims(Vols_tkk);
SizeData=size(Vols_tkk);
nVol=SizeData(end);

for vol=1:nVol
    Vols = fft(ifft(ifft(Vols_tkk,[],2),[],3),[],1);
end

Freq=MinMax_pt(1):MinMax_pt(2);
s=sprintf('%s_Comparison_LipRem.ps',name_fig);
delete(s); 
f=0;

vox_step=round(sqrt(sum(Mask2dBrain(:))/30));
index=0;
 for r=1:vox_step:SizeData(2);
        for c=1:vox_step:SizeData(3);
            if Mask2dBrain(r,c)
               index=index+1;
               SelVox(index,1)=r;
               SelVox(index,2)=c;
            end
        end
 end

cc=hsv(nVol);
ImSiR=20480;ImSiC=40960;
    for i=1:size(SelVox,1);
        f=f+1;
        figs(f)=figure('visible', 'off','Position', [100, 200, ImSiC,ImSiR]); 
        subplot(1,2,1),
        hold on
        for vol=1:nVol 
        	plot(Freq,abs(Vols(Freq,SelVox(i,1),SelVox(i,2),vol)),'color',cc(vol,:))
        end
        hold off
        h_legend=legend(legends,'Location','northwest','boxoff');
        set(h_legend,'FontSize',4);
        axis square;
       print(figs(f), '-append', '-dpsc2', s); 
        
        f=f+1;
        figs(f)=figure('visible', 'off','Position', [100, 200, ImSiC,ImSiR]); 
       % image2plot=squeeze(sum(abs(reshUres(:,:,2,k)),3));
       if nDims ==4
           image2plot=squeeze(sum(sum(abs(Vols(:,:,:,1)),1),4));
       elseif nDims ==3
           image2plot=squeeze(sum(abs(Vols(:,:,:,1)),1));
       end
        image2plot(SelVox(i,1),SelVox(i,2))=2*(max(image2plot(:)));
        
        subplot(1,2,2),imagesc(image2plot/(max(image2plot(:))));  
        axis('off');
        axis square;
        colormap default 
      
        print(figs(f), '-append', '-dpsc2', s); 

    end;
close all;




end

