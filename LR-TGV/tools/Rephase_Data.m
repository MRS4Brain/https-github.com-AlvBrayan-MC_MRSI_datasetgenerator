function Filtered_Data_trr=Rephase_Data(Data_trr,BrainMask2D)
Filtered_Data_trr=zeros(size(Data_trr));
 for a=1:size(Data_trr,2)
     for b=1:size(Data_trr,3)
         if(BrainMask2D(a,b))    
            %  Filtered_Data_trr(:,a,b)=KKRecurs(Data_trr(:,a,b)); 
               Filtered_Data_trr(:,a,b)=Compute_MinPhase(Data_trr(:,a,b)); 
         end
     end
 end
end
