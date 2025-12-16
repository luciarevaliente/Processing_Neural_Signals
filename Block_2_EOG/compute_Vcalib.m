function[Vcal_minus,Vcal_plus] = compute_Vcalib(EOG)
% compute_Vcalib
% This function determines Vcalib- and Vcalib+ (Task 6) from a calibration
% sequence (see lab manual appendix)

EOG_sign = sign(EOG);
EOG_sign = abs(EOG_sign(2:end)-EOG_sign(1:end-1));

zc_idx = find(EOG_sign);

Vtmp = zeros(length(zc_idx)-1,1);
for k = 1:length(Vtmp)
   Vtmp(k) = sign(EOG(round(mean([zc_idx(k),zc_idx(k+1)]))))...
       *abs(max(EOG(zc_idx(k):zc_idx(k+1)))); 
end

Vcal_minus = mean(Vtmp(Vtmp<0));
Vcal_plus = mean(Vtmp(Vtmp>0));