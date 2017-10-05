function [peakWeightDiff, peakWeightRatio] = cal_pTWpeakWeight(conA, conB)

for i = 1:length(conA)
    
   a(i) = conA(i,i);
   b(i) = conB(i,i);
   
   peakWeightDiff(i) = b(i) - a(i);
   
   peakWeightRatio(i) = b(i) / a(i);
    
end


end