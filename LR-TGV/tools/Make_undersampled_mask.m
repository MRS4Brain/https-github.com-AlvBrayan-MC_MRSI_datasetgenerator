function [US_MASK]=Make_undersampled_mask(expo,Fil_Fact,kmask,safe_radius)

% kmask=fftshift(fftshift(fftshift(kmask,1),2),3);
Size_grid=size(kmask);
for a=1:Size_grid(1);
    for b=1:Size_grid(2);
        %radius(a,b)=sqrt((a-Size_grid(1)*0.5-1)^2+ (b-Size_grid(2)*0.5-1)^2)+1; 
        radius(a,b)=(a-Size_grid(1)*0.5-1)^2+ (b-Size_grid(2)*0.5-1)^2 -safe_radius^2;
        if radius(a,b)<0;radius(a,b)=0;else
        radius(a,b)=sqrt( radius(a,b)/(2*Size_grid(1)*0.5*Size_grid(2)*0.5-safe_radius^2) );
        end
    end
end



%Prob_Map=((radius/safe_radius).^(expo));% way 1 to insure  safe_radius
%Prob_Map((end),round(end))

Prob_Map=(1-radius).^expo;

if expo==0
    Prob_Map=double(radius==0);
end

US_MASK=fftshift(fftshift(kmask,1),2);
% US_MASK=kmask;
Tot_AcqP=sum(kmask(:));
fil=1;
while(fil>Fil_Fact)
    a=round(rand()*(Size_grid(1)-1))+1;
    b=round(rand()*(Size_grid(2)-1))+1;
     if(Prob_Map(a,b)<rand()) ;
        US_MASK(a,b)=0;
     end
    fil=sum(US_MASK(:))/Tot_AcqP;
end

% US_MASK=fftshift(fftshift(fftshift(US_MASK,1),2),3);
US_MASK=fftshift(fftshift(US_MASK,1),2);
%imagesc(kmask)
figure
imagesc(US_MASK)
end