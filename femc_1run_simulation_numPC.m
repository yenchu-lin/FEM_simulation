% Test the weights of individual PCs
% Simulate LGN activity and behavior responses during letter viewing
% Last update on 2022/04/12 by YCL 

% Inputs
disp('Test the weights of individual PCs')

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
info.img.nimg = 100; % 100;

% 3. noise parameters
info.img.noisetype = '1/F'; % 1/F or Gaussian or none
info.img.cnoise = 0.3;
info.img.cstim  = 0.0002; % 0.0003
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
info.eyemvnt.type = 'exp_traj'; % 'sim_traj'; % 'exp_traj' or 'sim_traj'
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
info.disc.pc = 20; % number of principal components

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

img_size = info.img.imsize;
emv_x = info.emv.x(:,1,itraj,icell);
emv_y = info.emv.x(:,2,itraj,icell);
t = 0 : info.eyemvnt.conversionFactor : (length(emv_x)-1)*info.eyemvnt.conversionFactor;
ts_sfilt = squeeze(ts_spatial(:,1,:) - ts_spatial(:,2,:)); % center - surround

% PCA
% get the weight of each pc
[nimg,nt,nstim] = size(TS);
data = [];
for istim = 1:nstim
    data = [data;TS(:,:,istim)];
end

[coeff,score,latent,tsquared,explained,mu] = pca(data);

%% plot the weights of each pc
iplot = 15;
figure('Position',[10 10,800,600])
plot(explained(1:iplot),'o-','linewidth',2)
xlabel('principal components')
ylabel('variances')
title('variances of individual principal components')
set(gca,'fontsize',18,'linewidth',2)
box off

% plot - save figure
saveplot = 1;
fname = 'ch2_pc_weights';
figpt = '/Users/mb/Desktop/Thesis/Figures';
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end
