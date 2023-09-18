%WeightedAvg_SEM.m

function [wtdavg,wtdSEM,wtdSD] = WeightedAvg_SEM(data,wts)

wtdavg = sum(data.*wts)/sum(wts);

n = numel(data);
wtdx2 = sum(wts.*data.^2)/sum(wts);
wtd2mom = wtdx2 - wtdavg^2;
df_factor = (sum(wts)^2)/((sum(wts)^2) - sum(wts.^2));
df_a = (sum(wts)^2);
df_b = sum(wts.^2);
wtdVar = wtd2mom*df_factor;
wtdSD = sqrt(wtdVar);
wtdSEM = sqrt(wtd2mom*(df_b/(df_a-df_b)));

% k_stim = abs(CombinedKFmat(:,2));
% dk = CombinedKFmat(:,3); %SD k
% wts = 1./(dk).^2;
% speed_stim = 2*pi*0.1./k_stim;
% dv = speed_stim.*dk./k_stim; %SD c
% cwts = 1./(dv).^2; %1/variance
% n = numel(k_stim);
% 
% %|k| mean and SE (weighted)
% k_wtdmean = sum(k_stim.*wts)/sum(wts)
% wtdx2 = sum(wts.*k_stim.^2)/sum(wts);
% wtd2mom = wtdx2 - k_wtdmean^2;
% df_factor = (sum(wts)^2)/((sum(wts)^2) - sum(wts.^2));
% df_a = (sum(wts)^2);
% df_b = sum(wts.^2);
% kwtdVar = wtd2mom*df_factor;
% kwtdSEM = sqrt(wtd2mom*(df_b/(df_a-df_b)))

end