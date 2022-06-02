% Plotting firing patterns for thesis figures
% Simulate LGN activity and behavior responses during letter viewing
% Last updated on 2022/04/12

% 1. Simulated eye movements (move horizontally from left to right) for
% letter viewing (Figure 5)
% 2. Measured human eye movements for letter viewing - plots multiple
% firings (Figure 6)
% 3. Measured human eye movements for letter viewing without noise (Figure
% 6)

%% 1. Simulated eye movements (move horizontally from left to right) for letter viewing 
% with and without noises but without ramp
saveplot = 0;
disp('simulated eye movements and letters - for thesis figures')

% 1. input stimuli 
% 1.1. input type - letters
info.img.type = 'letter';
info.img.nstim = 2;
info.stim = cell(info.img.nstim,1);
info.stim{1} = 'H';
info.stim{2} = 'N';

% 2.image parameters
info.img.imsize = 512;
info.img.D = 2;
info.img.npad = 2;
info.img.nimg = 1; % 100;

% 3. noise parameters
info.img.noisetype = 'none'; % 1/F or Gaussian or none
info.img.cnoise = 0.3;
info.img.cstim  = 0.0005; % 0.0003; % 0.0003
% 3.1 Gaussian Noise
info.img.varnoise = 1;
info.img.cutoffnoise = 2;
info.img.shift = 0;

% 3.2 1/F noise parameters
info.img.freq = 20; % for 1/F noise mask 
info.img.noisefreqlimit = 2^6; % upper bound frequency limit of the 1/f noise (range: 1-noisefreqlimit)
info.img.noisefrequp = 2^9;
info.img.noisefreqlo = 2^0;

% 4. eye movement parameters
info.eyemvnt.ddpi = 1;
info.eyemvnt.type = 'sim_traj'; %'sim_traj'; % 'exp_traj' or 'sim_traj'
info.eyemvnt.datadir = '/Users/mb/Documents/Data/Femc_eyetraj/'; % '/Users/mb/Documents/Data/Femc_simulation/eye_traj';
info.eyemvnt.filename = 'Data_JI.mat'; % info.eyemvnt.filename = 'clean_Christie_8cd.mat';
info.eyemvnt.fname = 'c1'; % c1 (HN) or c2 (EF) based on the dataset
info.eyemvnt.rotate = 0; % in degrees to rotate the eye trajectories
info.eyemvnt.t.baseline = 0; % time information of the ramp protocol
info.eyemvnt.t.ramp = 1000;    % time information of the ramp protocol
info.eyemvnt.t.plateau = 500;  % time information of the ramp protocol
info.eyemvnt.t.post = 0;       % time information of the ramp protocol
% 
% % 5. neural population parameters
info.cell.ncell = 5; % number of cells
info.cell.type = 'P'; % P or M cell
info.cell.rftype = 'ON'; % ON or OFF cell
info.cell.spat.model_no = 1; % 1- circular Gaussian
info.cell.temp.delta_t = 0.001;
info.cell.temp.model_no = 3; % 2- gamma; 3- HPLP

% 6. discriminator parameters
info.disc.train_percent = 0.80; % percentage of data into the training 
info.disc.pc = 2; % number of principal components

% Simulate LGN activity during letter viewing (2019.08.19)
[model,info] = femc_lgn_parameters(info);

% randomized cell center - use "femc_cell_center"
% [info] = femc_cell_center(info,model);
% for testing purpose; can only have up to 63 cells
[info] = femc_cell_center_old(info,model); 

[info] = femc_input_eye_movement(info);

% istim/itraj/icell/iimg
imsize = info.img.imsize * info.img.npad;
fimg = zeros(imsize,imsize,info.img.nstim,2,info.img.nimg);
for iimg = 1:info.img.nimg
    [img,img_rotate,noiseMask] = femc_stimulus(info);
    [filterimg] = femc_spatial_filter_img(info,img,model);
    fimg(:,:,:,:,iimg) = filterimg;
end

icell = 1;
itraj = 1;
for iimg = 1:info.img.nimg
    [ts_spatial] = femc_spatial(info,fimg(:,:,:,:,iimg),info.emv.x(:,:,itraj,icell)); 
%             [ts_ramp,ramp_contrast] = femc_ramp(ts_spatial,info);
    [Timeseries,info] = femc_temporal(model,ts_spatial,info);
%             [Timeseries,info] = femc_temporal(model,ts_ramp,info);
    TS(iimg,:,:) = Timeseries;
end
disp('traj')
disp(itraj)
disp('cell')
disp(icell)

img_size = info.img.imsize;
emv_x = info.emv.x(:,1,itraj,icell);
emv_y = info.emv.x(:,2,itraj,icell);
t = 0 : info.eyemvnt.conversionFactor : (length(emv_x)-1)*info.eyemvnt.conversionFactor;
ts_sfilt = squeeze(ts_spatial(:,1,:) - ts_spatial(:,2,:)); % center - surround

%% plot single condition
% plot the eye trajectories on the letters & firing rates
ylimmax = [5*10^-4, 0.02];
femc_plot_letter_firing(ts_sfilt,TS,img,info,emv_x,emv_y,t,ylimmax)

% save figure
if strcmp(info.img.noisetype,'1/F')
    fname = ['ch2_letter_firing_' info.cell.type '_cell_noise_1F'];
else
    fname = ['ch2_letter_firing_' info.cell.type '_cell_noise_' info.img.noisetype];
end
figpt = '/Users/mb/Desktop/Thesis/Figures';
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end

%% 2. Measured human eye movements for letter viewing - plots multiple firings
% with noises, with ramp
saveplot = 0;
disp('human eye movements and letters - for thesis figures')

% 1. input stimuli 
% 1.1. input type - letters
info.img.type = 'letter';
info.img.nstim = 2;
info.stim = cell(info.img.nstim,1);
info.stim{1} = 'H';
info.stim{2} = 'N';

% 2.image parameters
info.img.imsize = 512;
info.img.D = 2;
info.img.npad = 2;
info.img.nimg = 10; % 100;

% 3. noise parameters
info.img.noisetype = '1/F'; % 1/F or Gaussian or none
info.img.cnoise = 0.3;
info.img.cstim  = 0.0003; % 0.0003; % 0.0003
% 3.1 Gaussian Noise
info.img.varnoise = 1;
info.img.cutoffnoise = 2;
info.img.shift = 0;

% 3.2 1/F noise parameters
info.img.freq = 20; % for 1/F noise mask 
info.img.noisefreqlimit = 2^6; % upper bound frequency limit of the 1/f noise (range: 1-noisefreqlimit)
info.img.noisefrequp = 2^9;
info.img.noisefreqlo = 2^0;

% 4. eye movement parameters
info.eyemvnt.ddpi = 1;
info.eyemvnt.type = 'exp_traj'; %'exp_traj'; % 'exp_traj' or 'sim_traj'
info.eyemvnt.datadir = '/Users/mb/Documents/Data/Femc_eyetraj/'; % '/Users/mb/Documents/Data/Femc_simulation/eye_traj';
info.eyemvnt.filename = 'Data_JI.mat'; % info.eyemvnt.filename = 'clean_Christie_8cd.mat';
info.eyemvnt.fname = 'c1'; % c1 (HN) or c2 (EF) based on the dataset
info.eyemvnt.rotate = 0; % in degrees to rotate the eye trajectories
info.eyemvnt.t.baseline = 0; % time information of the ramp protocol
info.eyemvnt.t.ramp = 1000;    % time information of the ramp protocol
info.eyemvnt.t.plateau = 500;  % time information of the ramp protocol
info.eyemvnt.t.post = 0;       % time information of the ramp protocol
% 
% % 5. neural population parameters
info.cell.ncell = 5; % number of cells
info.cell.type = 'P'; % P or M cell
info.cell.rftype = 'ON'; % ON or OFF cell
info.cell.spat.model_no = 1; % 1- circular Gaussian
info.cell.temp.delta_t = 0.001;
info.cell.temp.model_no = 3; % 2- gamma; 3- HPLP

% 6. discriminator parameters
info.disc.train_percent = 0.80; % percentage of data into the training 
info.disc.pc = 2; % number of principal components

% Simulate LGN activity during letter viewing (2019.08.19)
[model,info] = femc_lgn_parameters(info);

% randomized cell center - use "femc_cell_center"
% [info] = femc_cell_center(info,model);
% for testing purpose; can only have up to 63 cells
[info] = femc_cell_center_old(info,model); 

[info] = femc_input_eye_movement(info);

% istim/itraj/icell/iimg
imsize = info.img.imsize * info.img.npad;
fimg = zeros(imsize,imsize,info.img.nstim,2,info.img.nimg);
for iimg = 1:info.img.nimg
    [img,img_rotate,noiseMask] = femc_stimulus(info);
    [filterimg] = femc_spatial_filter_img(info,img,model);
    fimg(:,:,:,:,iimg) = filterimg;
end

TS = zeros(info.img.nimg,info.eyemvnt.nsamp,info.img.nstim);
icell = 1;
itraj = 1;
for iimg = 1:info.img.nimg
    [ts_spatial] = femc_spatial(info,fimg(:,:,:,:,iimg),info.emv.x(:,:,itraj,icell)); 
    [ts_ramp,ramp_contrast] = femc_ramp(ts_spatial,info);
    [Timeseries,info] = femc_temporal(model,ts_ramp,info);
    TS(iimg,:,:) = Timeseries;
end
disp('traj')
disp(itraj)
disp('cell')
disp(icell)

img_size = info.img.imsize;
emv_x = info.emv.x(:,1,itraj,icell);
emv_y = info.emv.x(:,2,itraj,icell);
t = 0 : info.eyemvnt.conversionFactor : (length(emv_x)-1)*info.eyemvnt.conversionFactor;

%% plot1 - firing from multiple trajs
% istim = 1;
% itraj = 2;
nsamp = info.eyemvnt.nsamp;
ntraj = info.eyemvnt.ntraj;
nimg = info.img.nimg;
% figure; plot(TS(iimg,:,istim,itraj))
TS1 = squeeze(TS(:,:,1));
TS2 = squeeze(TS(:,:,2));
 
figure
plot3(repmat(t,[nimg,1])', repmat(1:nimg,[nsamp,1]),TS1, 'linewidth',1.5);
set(gca,'fontsize',18,'linewidth',2)
zlim([0 0.015])
xlabel('time (ms)')
ylabel('images')
zlabel('firing rate')

% plot - save figure
fname = ['ch2_letter_firing_multicells' info.cell.type];
figpt = '/Users/mb/Desktop/Thesis/Figures';
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end

%% plot2 - single condition
% plot the eye trajectories on the letters & firing rates
ylimmax = [4*10^-4, 0.01];
ts_sfilt = squeeze(ts_ramp(:,1,:) - ts_ramp(:,2,:)); % center - surround
femc_plot_letter_firing(ts_sfilt,TS,img,info,emv_x,emv_y,t,ylimmax)

% save figure
if strcmp(info.img.noisetype,'1/F')
    fname = ['ch2_letter_firing_' info.cell.type '_cell_noise_1F_humanEM'];
else
    fname = ['ch2_letter_firing_' info.cell.type '_cell_noise_' info.img.noisetype '_humanEM'];
end
figpt = '/Users/mb/Desktop/Thesis/Figures/ch2';
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end
 
%% Plot3 - one specfic traj on the image
istim = 1;
figure('Position', [10 10 1200 500])
subplot(1,2,1)
imagesc(flipud(img(:,:,istim)))
set (gca,'Ydir','normal')
colormap(gray)
hold on
t = 0:2.9070:1501;

plot(emv_x,emv_y,'r-','linewidth',2) %,'color',cmap)
axis square
ticksize = info.img.D*2;
imsize = info.img.imsize*info.img.npad; % in pixel

% label in degrees
A=(-(info.img.D*info.img.npad/2) : info.img.D*info.img.npad/ticksize : info.img.D*info.img.npad/2)';
s = num2str(A);
set(gca,'XTick',(0:1/ticksize:1)*imsize,'XTickLabel',s)
set(gca,'YTick',(0:1/ticksize:1)*imsize,'YTickLabel',s)
xlabel('x (deg)')
ylabel('y (deg)') 
set(gca, 'FontSize', 16)

subplot(1,2,2)
plot(t,emv_x-imsize/2,'linewidth',2)
hold on 
plot(t,emv_y-imsize/2,'linewidth',2)
legend('x','y','location','southeast')
xlabel('time (ms)')
ylabel('position (arcmin)') 
set(gca,'FontSize', 16)
set(gca,'linewidth',2)
box off

% plot - save figure
fname = ['ch2_letter_firing_eye_traj_' info.cell.type '_humanEM'];
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end

%% 3. Measured human eye movements for letter viewing - plots multiple firings - no noise
% with noises, with ramp
saveplot = 0;
disp('human eye movements and letters - no noise - for thesis fig 6 - no noise')

% 1. input stimuli 
% 1.1. input type - letters
info.img.type = 'letter';
info.img.nstim = 2;
info.stim = cell(info.img.nstim,1);
info.stim{1} = 'H';
info.stim{2} = 'N';

% 2.image parameters
info.img.imsize = 512;
info.img.D = 2;
info.img.npad = 2;
info.img.nimg = 10; % 100;

% 3. noise parameters
info.img.noisetype = 'none'; % 1/F or Gaussian or none
info.img.cnoise = 0.3;
info.img.cstim  = 0.0003; % 0.0003; % 0.0003
% 3.1 Gaussian Noise
info.img.varnoise = 1;
info.img.cutoffnoise = 2;
info.img.shift = 0;

% 3.2 1/F noise parameters
info.img.freq = 20; % for 1/F noise mask 
info.img.noisefreqlimit = 2^6; % upper bound frequency limit of the 1/f noise (range: 1-noisefreqlimit)
info.img.noisefrequp = 2^9;
info.img.noisefreqlo = 2^0;

% 4. eye movement parameters
info.eyemvnt.ddpi = 1;
info.eyemvnt.type = 'exp_traj'; %'exp_traj'; % 'exp_traj' or 'sim_traj'
info.eyemvnt.datadir = '/Users/mb/Documents/Data/Femc_eyetraj/'; % '/Users/mb/Documents/Data/Femc_simulation/eye_traj';
info.eyemvnt.filename = 'Data_JI.mat'; % info.eyemvnt.filename = 'clean_Christie_8cd.mat';
info.eyemvnt.fname = 'c1'; % c1 (HN) or c2 (EF) based on the dataset
info.eyemvnt.rotate = 0; % in degrees to rotate the eye trajectories
info.eyemvnt.t.baseline = 0; % time information of the ramp protocol
info.eyemvnt.t.ramp = 1000;    % time information of the ramp protocol
info.eyemvnt.t.plateau = 500;  % time information of the ramp protocol
info.eyemvnt.t.post = 0;       % time information of the ramp protocol
% 
% % 5. neural population parameters
info.cell.ncell = 5; % number of cells
info.cell.type = 'P'; % P or M cell
info.cell.rftype = 'ON'; % ON or OFF cell
info.cell.spat.model_no = 1; % 1- circular Gaussian
info.cell.temp.delta_t = 0.001;
info.cell.temp.model_no = 3; % 2- gamma; 3- HPLP

% 6. discriminator parameters
info.disc.train_percent = 0.80; % percentage of data into the training 
info.disc.pc = 2; % number of principal components

% Simulate LGN activity during letter viewing (2019.08.19)
[model,info] = femc_lgn_parameters(info);

% randomized cell center - use "femc_cell_center"
% [info] = femc_cell_center(info,model);
% for testing purpose; can only have up to 63 cells
[info] = femc_cell_center_old(info,model); 

[info] = femc_input_eye_movement(info);

% istim/itraj/icell/iimg
imsize = info.img.imsize * info.img.npad;
fimgnn = zeros(imsize,imsize,info.img.nstim,2,info.img.nimg);
for iimg = 1:info.img.nimg
    [imgnn,img_rotatenn,noiseMasknn] = femc_stimulus(info);
    [filterimgnn] = femc_spatial_filter_img(info,img,model);
    fimgnn(:,:,:,:,iimg) = filterimgnn;
end

TS1_nonoise = zeros(info.img.nimg,info.eyemvnt.nsamp,info.img.nstim);
icell = 1;
itraj = 1;
for iimg = 1:info.img.nimg
    [ts_spatial] = femc_spatial(info,fimgnn(:,:,:,:,iimg),info.emv.x(:,:,itraj,icell)); 
    [ts_ramp,ramp_contrast] = femc_ramp(ts_spatial,info);
    [Timeseries,info] = femc_temporal(model,ts_ramp,info);
    TS1_nonoise(iimg,:,:) = Timeseries;
end
disp('traj')
disp(itraj)
disp('cell')
disp(icell)

img_size = info.img.imsize;
emv_x = info.emv.x(:,1,itraj,icell);
emv_y = info.emv.x(:,2,itraj,icell);
t = 0 : info.eyemvnt.conversionFactor : (length(emv_x)-1)*info.eyemvnt.conversionFactor;

%% plot1 - firing from multiple trajs
nsamp = info.eyemvnt.nsamp;
ntraj = info.eyemvnt.ntraj;
nimg = info.img.nimg;
nimg = 1;
% figure; plot(TS(iimg,:,istim,itraj))
TS1 = squeeze(TS1_nonoise(:,:,1));
TS2 = squeeze(TS1_nonoise(:,:,2));
 
figure
plot3(repmat(t,[nimg,1])', repmat(1:nimg,[nsamp,1]),TS1, 'linewidth',1.5);
set(gca,'fontsize',18,'linewidth',2)
zlim([0 0.015])
xlabel('time (ms)')
ylabel('images')
zlabel('firing rate')

% plot - save figure
fname = ['ch2_letter_firing_multicells' info.cell.type 'humanEM_nonoise'];
figpt = '/Users/mb/Desktop/Thesis/Figures/ch2';
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end

% plot noise and no noise together
% nsamp = info.eyemvnt.nsamp;
% ntraj = info.eyemvnt.ntraj;
% nimg = info.img.nimg;
% nimg = nimg+1;
% % figure; plot(TS(iimg,:,istim,itraj))
% TS_TS_nonoise = [TS1_nonoise(1,:,:);TS];
% TS1 = squeeze(TS_TS_nonoise(:,:,1));
% % TS2 = squeeze(TS(:,:,1));
%  
% figure
% plot3(repmat(t,[nimg,1])', repmat(1:nimg,[nsamp,1]),TS1, 'linewidth',1.5);
% set(gca,'fontsize',18,'linewidth',2)
% xlabel('time (ms)')
% ylabel('images')
% zlabel('firing rate')
% 
% % plot - save figure
% fname = ['ch2_letter_firing_multicells' info.cell.type 'humanEM_bothnoise'];
% figpt = '/Users/mb/Desktop/Thesis/Figures/ch2';
% if saveplot
%     saveas(gca,fullfile(figpt,fname), 'tif')
% end

%% Plot - one specfic traj on the image
istim = 1;
figure('Position', [10 10 500 500])
imagesc(flipud(img_nonoise(:,:,istim)))
set (gca,'Ydir','normal')
colormap(gray)
caxis([-0.0012, 0.0013]);

hold on
t = 0:2.9070:1501;

plot(emv_x,emv_y,'r-','linewidth',2) %,'color',cmap)
axis square
ticksize = info.img.D*2;
imsize = info.img.imsize*info.img.npad; % in pixel

% label in degrees
A=(-(info.img.D*info.img.npad/2) : info.img.D*info.img.npad/ticksize : info.img.D*info.img.npad/2)';
s = num2str(A);
set(gca,'XTick',(0:1/ticksize:1)*imsize,'XTickLabel',s)
set(gca,'YTick',(0:1/ticksize:1)*imsize,'YTickLabel',s)
xlabel('x (deg)')
ylabel('y (deg)') 
set(gca, 'FontSize', 16)

fname = ['ch2_letter_firing_' info.cell.type '_cell_noise_' info.img.noisetype '_humanEM_nonoise'];
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end

%% plot2 - single condition
% plot the eye trajectories on the letters & firing rates
ylimmax = [2*10^-4, 0.01];
ts_sfilt = squeeze(ts_ramp(:,1,:) - ts_ramp(:,2,:)); % center - surround
femc_plot_letter_firing(ts_sfilt,TS,img,info,emv_x,emv_y,t,ylimmax)

% save figure
if strcmp(info.img.noisetype,'1/F')
    fname = ['ch2_letter_firing_' info.cell.type '_cell_noise_1F_humanEM_nonoise'];
else
    fname = ['ch2_letter_firing_' info.cell.type '_cell_noise_' info.img.noisetype '_humanEM_nonoise'];
end
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end

%% Plot3 - one specfic traj on the image
istim = 1;
figure('Position', [10 10 1200 500])
subplot(1,2,1)
imagesc(flipud(img(:,:,istim)))
set (gca,'Ydir','normal')
colormap(gray)
hold on
t = 0:2.9070:1501;

plot(emv_x,emv_y,'r-','linewidth',2) %,'color',cmap)
axis square
ticksize = info.img.D*2;
imsize = info.img.imsize*info.img.npad; % in pixel

% label in degrees
A=(-(info.img.D*info.img.npad/2) : info.img.D*info.img.npad/ticksize : info.img.D*info.img.npad/2)';
s = num2str(A);
set(gca,'XTick',(0:1/ticksize:1)*imsize,'XTickLabel',s)
set(gca,'YTick',(0:1/ticksize:1)*imsize,'YTickLabel',s)
xlabel('x (deg)')
ylabel('y (deg)') 
set(gca, 'FontSize', 16)

subplot(1,2,2)
plot(t,emv_x-imsize/2,'linewidth',2)
hold on 
plot(t,emv_y-imsize/2,'linewidth',2)
legend('x','y','location','southeast')
xlabel('time (ms)')
ylabel('position (arcmin)') 
set(gca,'FontSize', 16)
set(gca,'linewidth',2)
box off

% plot - save figure
fname = ['ch2_letter_firing_eye_traj_' info.cell.type '_humanEM_nonoise'];
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end
