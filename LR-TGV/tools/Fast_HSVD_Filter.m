function [filtMRSI ] = Fast_HSVD_Filter(mrsiData,tempSamplerate,modelOrder,Water_minfreq,Water_maxfreq)

minFreq         = Water_minfreq;%in hz
maxFreq         = Water_maxfreq;%in hz
filtMRSI = zeros(size(mrsiData));

[freqs,~,basis,amps] = HSVD(mrsiData, tempSamplerate, modelOrder);
indx                 = find((freqs >= minFreq) & (freqs <= maxFreq));
filtMRSI = mrsiData-sum(basis(:, indx) * diag(amps(indx), 0), 2);

end

