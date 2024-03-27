function out = CombEncodOp(Data_rr,kmask,SENSE)

out= AdjEncodOp(ForwEncodOp(Data_rr,kmask,SENSE), kmask,SENSE);

end