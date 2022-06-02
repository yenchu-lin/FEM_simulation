function [] = femc_plot_grat_sf_tuning(auc_spat,auc_temp,sf_list,info,sscale)
% plot sensitivity v1 - by area under the curve

figure('Position',[10 10 600 600])
hold on; 
plot(sf_list,auc_spat(1,:)*sscale,'o-','linewidth',2)
plot(sf_list,auc_temp(1,:),'o-','linewidth',2)
xlabel('spatial freq (cpd)')
ylabel('area under the firing rate curve')
legend('spatial','spatial+temporal')
title([info.cell.type, ' cell'])
set(gca,'Fontsize',18,'linewidth',2)

end

