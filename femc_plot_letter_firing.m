function [] = femc_plot_letter_firing(ts_sfilt,TS,img,info,emv_x,emv_y,t,ylimmax)
% plot single condition
% plot the eye trajectories on the letters & firing rates

iimg = 1;
ts_sfilt(ts_sfilt<0) = 0;

figure('Position',[2300,10,800,1400])
for iletter = 1:2
    
%     subplot(2,3,1+(iletter-1)*3)
    subplot(3,2,iletter)
    hold on
    imagesc(img(:,:,iletter))
    colormap(gray)
    plot(emv_x, emv_y,'linewidth',1.5)
%     img_size = info.img.imsize;
%     xlim([img_size/2+1 img_size+img_size/2])
%     ylim([img_size/2+1 img_size+img_size/2])
    xlim([0 info.img.imsize*info.img.npad])
    ylim([0 info.img.imsize*info.img.npad])
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    axis square
    set(gca, 'YDir','reverse')
    
%     subplot(2,3,2+(iletter-1)*3)
    subplot(3,2,2+iletter)
    plot(t,ts_sfilt(:,iletter),'linewidth',1.5)
    ylim([0 ylimmax(1)])
    title('spatial filter')
    if iletter == 1
        ylabel('firing profile')
    end
    if iletter == 2
        set(gca,'yticklabel',{[]})
    end
    set(gca,'linewidth',2,'fontsize',18)
    box off
    
%     subplot(2,3,3+(iletter-1)*3)
    subplot(3,2,4+iletter)
    hold on
    plot(t,TS(iimg,:,iletter),'linewidth',1.5)
    ylim([0 ylimmax(2)])
    xlabel('time(ms)')
    title('spatiotemproal filter')
    if iletter == 1
        ylabel('firing profile')
    end   
    if iletter == 2
        set(gca,'yticklabel',{[]})
    end
    set(gca,'linewidth',2,'fontsize',18)
    box off
end
end

