function [] = femc_plot_grat_sf_firing(grat_img,fr_spat,fr_temp,emv_x,emv_y,sf_list,info,t,nsf,img_size,ymax)
% plot sf vs firing 

iletter = 1;
figure('Position',[10 10 2600 800])
for itest = 1:nsf
    subplot(3,nsf,itest)
    imagesc(grat_img(:,:,iletter,itest))
    colormap(gray)
    xlim([img_size/2+1 img_size+img_size/2])
    ylim([img_size/2+1 img_size+img_size/2])
    axis square
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    title([num2str(sf_list(itest)), ' cpd'])
    if itest == 1
        hold on
        plot(emv_x, emv_y)
        ylabel([info.cell.type, ' cell'])
    end
    
    subplot(3,nsf,itest+nsf)
    plot(t,fr_spat(:,iletter,itest),'linewidth',2)
    ylim([0 ymax(1)])
    if itest == 1
        ylabel('spatial filter')
    else
        set(gca,'xticklabel',{[]})
        set(gca,'yticklabel',{[]})
    end
    set(gca,'Fontsize',18,'linewidth',2)
    box off
    
    subplot(3,nsf,itest+2*nsf)
    plot(t,fr_temp(:,iletter,itest),'linewidth',2)
    ylim([0 ymax(2)])
    if itest == 1
        ylabel('spatial+temporal filter')
    else
        set(gca,'xticklabel',{[]})
        set(gca,'yticklabel',{[]})
    end
    set(gca,'Fontsize',18,'linewidth',2)
    box off
end
end

