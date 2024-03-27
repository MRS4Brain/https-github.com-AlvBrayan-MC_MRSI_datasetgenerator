function encodes=MakeEllipEncoding(d1,d2,d3)
    
    encodes=zeros(d1,d2,d3);
    
    
    lowd1 = -floor(d1/2);
      if(mod(d1,2))  uppd1 =floor(d1/2);
      else uppd1 =floor(d1/2)-1;
      end;
     lowd2 = -floor(d2/2);
      if(mod(d2,2))  uppd2 =floor(d2/2);
      else uppd2 =floor(d2/2)-1;
      end;
     lowd3 = -floor(d3/2);
     if(mod(d3,2)) uppd3 = floor(d3/2);
     else uppd3 = floor(d3/2)-1;
     end;

	
    addi_voxel_rad=0.125;
     for a3=lowd3:uppd3
         for a2=lowd2:uppd2
             for a1=lowd1:uppd1
        
	                   
						  if(uppd1 == 0)  d =0;
                          else d =(a1/(uppd1+addi_voxel_rad));
                          end;
						 dist = d*d;
						 if(uppd2 == 0)  d=0;
                         else d =(a2/(uppd2+addi_voxel_rad));
                         end;
						 dist = dist+ d*d;
						 if(uppd3 == 0)  d=0;
                         else d =(a3/(uppd3+addi_voxel_rad));
                         end;
						 dist = dist+d*d;
						 dist = sqrt( dist );
	                     
						 if (dist <= 1 )
                            encodes(a1-lowd1+1,a2-lowd2+1,a3-lowd3+1) =1; 
                         end
                         
         end;end;end;
   