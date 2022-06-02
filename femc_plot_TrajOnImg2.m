function [] = femc_plot_TrajOnImg2(data,icell,istim,itraj,cmap,isf,F,dn)
% The function plot one single trajectory on top of the letter, and its x
% and y against time
% for DDPI

figure('Position', [10 10 900 300])
subplot(1,2,1)
imagesc(flipud(data.img(:,:,istim)))
set (gca,'Ydir','normal')
colormap(gray)
hold on
x = data.info.emv.x(:,1,itraj,icell);
y = data.info.emv.x(:,2,itraj,icell);
t = 0:2.9070:1501;

% smooth parameter
% [sx,~] = sgfilt(x,isf,F,dn);
% [sy,~] = sgfilt(y,isf,F,dn);

plot(x,y,'r-','linewidth',2,'color',cmap)
axis square
ticksize = data.info.img.D*2;
imsize = data.info.img.imsize*data.info.img.npad; % in pixel

% label in degrees
A=(-(data.info.img.D*data.info.img.npad/2) : data.info.img.D*data.info.img.npad/ticksize : data.info.img.D*data.info.img.npad/2)';
s = num2str(A);
set(gca,'XTick',(0:1/ticksize:1)*imsize,'XTickLabel',s)
set(gca,'YTick',(0:1/ticksize:1)*imsize,'YTickLabel',s)
xlabel('x (deg)')
ylabel('y (deg)') 
set(gca, 'FontSize', 16)

subplot(1,2,2)
plot(t,x-imsize/2,'linewidth',2)
hold on 
plot(t,y-imsize/2,'linewidth',2)
legend('x','y')
xlabel('time (ms)')
ylabel('position (arcmin)') 
set(gca,'FontSize', 16)
set(gca,'linewidth',2)

end

