function [roi_pTW_av roi_pTW_recentredAv] = cal_ROI_pTW_av(betas_A,betas_B)

roi_pTW_av = mean([mean(betas_A,2),mean(betas_B,2)],2);

roi_pTW_recentredAv = [];

end