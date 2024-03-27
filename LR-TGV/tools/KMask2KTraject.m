function  kSpaceTrajectory=KMask2KTraject(kmask, FOV)
kdims=size(kmask);

is3D=(numel(size(kmask))>2);


pt=1;
if is3D
    for a=1:kdims(1)
        for b=1:kdims(2)
            for c=1:kdims(3)
                if kmask(a,b,c)
                    kSpaceTrajectory(1,pt)=(a-ceil((kdims(1)+1)/2))*2*pi/FOV(1);
                    kSpaceTrajectory(2,pt)=(b-ceil((kdims(2)+1)/2))*2*pi/FOV(2);
                    kSpaceTrajectory(3,pt)=(c-ceil((kdims(3)+1)/2))*2*pi/FOV(3);    
                pt=pt+1;
                end           

            end
        end
    end   
else
   for a=1:kdims(1)
        for b=1:kdims(2)
            if kmask(a,b)
                kSpaceTrajectory(1,pt)=(a-ceil((kdims(1)+1)/2))*2*pi/FOV(1);
                kSpaceTrajectory(2,pt)=(b-ceil((kdims(2)+1)/2))*2*pi/FOV(2);  
                pt=pt+1;
            end           

        end
    end  
end
