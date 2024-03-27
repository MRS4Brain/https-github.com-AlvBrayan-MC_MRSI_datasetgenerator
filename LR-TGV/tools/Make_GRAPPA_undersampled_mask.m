function [US_MASK]=Make_GRAPPA_undersampled_mask(Type,Fil_Fact,kmask)
kmask=fftshift(fftshift(kmask,1),2);

%Type imust be either "GRAPPA",  "CAIPI", "AAA"
Size_grid=size(kmask);
for a=1:Size_grid(1);
    for b=1:Size_grid(2);
        radius(a,b)=sqrt((a-Size_grid(1)*0.5-1)^2+ (b-Size_grid(2)*0.5-1)^2)+1;
    end
end

switch Type
			case 'GRAPPA'	
                grid=zeros(Size_grid);
                grid(1:2:end,1:2:end)=1;
                US_MASK=grid.*kmask;
            
            case 'CAIPI'
                
                grid=zeros(Size_grid);
                grid(1:4:end,1:2:end)=1;
                grid(3:4:end,2:2:end)=1;

                US_MASK=grid.*kmask;
            case 'AAA'
                
                grid=zeros(Size_grid);
                grid(1:4:end,1:4:end)=1;
                grid(2:4:end,3:4:end)=1;
                grid(3:4:end,2:4:end)=1;
                grid(4:4:end,4:4:end)=1;
                US_MASK=grid.*kmask;
                
end

                
Tot_AcqP=sum(kmask(:));
fil=0;
limit=1;
while(fil<Fil_Fact)
    US_MASK=US_MASK+(radius<limit);
    US_MASK(US_MASK>0)=1;
    fil=sum(US_MASK(:))/sum(kmask(:));
    limit=limit+0.1;
end
%{
figure
imagesc(kmask)
figure
imagesc(US_MASK)
%}
US_MASK=fftshift(fftshift(US_MASK,1),2);

end
