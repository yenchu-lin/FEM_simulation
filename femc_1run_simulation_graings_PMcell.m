%% Test gratings with different spatial frequencies
% Simulate LGN activity and behavior responses during letter viewing
% Last update on 2022/04/04 by YCL 

%% General inputs
% 1. input stimuli 
% input type - gratings 
info.img.type = 'grating';
info.img.nstim = 2;
% info.img.stimsf = sf_list(itest); % 5;
info.stim = cell(info.img.nstim,1);
info.stim{1} = 'vertical';
info.stim{2} = 'horizontal';
info.img.shift = 0; % shift the amount of the spatial period of the stimulus; one-quarter (0.25) or one-half (0.5)

% 2.image parameters
info.img.imsize = 512;
info.img.D = 2;
info.img.npad = 2;
info.img.nimg = 1; % 100;

% 3. noise parameters
info.img.noisetype = 'none'; % 1/F or Gaussian or none
info.img.cnoise = 0.3;
info.img.cstim  = 0.0003; % 0.0003
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
info.eyemvnt.type = 'sim_traj'; %'exp_traj'; % 'exp_traj' or 'sim_traj'
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
info.cell.rftype = 'ON'; % ON or OFF cell
info.cell.spat.model_no = 1; % 1- circular Gaussian
info.cell.temp.delta_t = 0.001;
info.cell.temp.model_no = 3; % 2- gamma; 3- HPLP
% info.cell.temp.celltype_no = 5; % 1- P ON cell or 3- P OFF cell or 5- M ON cell or 6 - M cell OFF cell
% info.cell.temp.model_no = 2;  % 2- gamma; 3- HPLP
% info.cell.temp.celltype_no = 1; % 1- P cell; 2- M cell

% 6. discriminator parameters
info.disc.train_percent = 0.80; % percentage of data into the training 
info.disc.pc = 2; % number of principal components
% info.disc.sbset.nocell = info.cell.ncell; % number of cells in each subset
% info.disc.sbset.no = 1; % randomly choose one index for the subset table

ntpoint = 1500;
sf_list = [.1,.2,.4,.8,1,2,3,4,5,6,7,8,9,10,12,16,24,32];
nsf = length(sf_list);

%% M cells simulation
info.cell.type = 'M'; % P or M cell
fr_spat = zeros(ntpoint,2,nsf);
fr_temp = zeros(ntpoint,2,nsf); 
auc_spat = zeros(2,nsf);
auc_temp = zeros(2,nsf);
max_spat = zeros(2,nsf);
max_temp = zeros(2,nsf);

icell = 1;
itraj = 1;
    
for itest = 1:nsf
    info.img.stimsf = sf_list(itest);
    
    % Simulate LGN activity during letter viewing (2019.08.19)
    % This simulation contains 5 parts
    %  1. Stimulus by femc_stimulus.m
    %  2. Spatial Filtering
    %  3. Eye movement (imcomplete: need to work on 'simulated_BW' option (2018.09.25))
    %  4. Ramp protocol for contrast
    %  5. Temporal Filtering
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
    grat_img(:,:,:,itest) = img;

    for iimg = 1:info.img.nimg
        [ts_spatial] = femc_spatial(info,fimg(:,:,:,:,iimg),info.emv.x(:,:,itraj,icell)); 
        [Timeseries,info] = femc_temporal(model,ts_spatial,info);
%       ignore ramp for grating simulation
%       [ts_ramp,ramp_contrast] = femc_ramp(ts_spatial,info);
%       [Timeseries,info] = femc_temporal(model,ts_ramp,info);
    end
    
    ts_sfilt = squeeze(ts_spatial(:,1,:) - ts_spatial(:,2,:)); % center - surround
    ts_sfilt(ts_sfilt<0) = 0;
    
    % save area under the curve 
    fr_spat(:,:,itest) = ts_sfilt;
    fr_temp(:,:,itest) = Timeseries;
    auc_spat(:,itest) = trapz(ts_sfilt);
    auc_temp(:,itest) = trapz(Timeseries);
    max_spat(:,itest) = max(ts_sfilt);
    max_temp(:,itest) = max(Timeseries);
    
    disp([num2str(sf_list(itest)) 'cpd done'])
end

%% plot 1 - sf vs firing 
ymax = [2*10^-4,0.1];
emv_x = info.emv.x(:,1,itraj,icell);
emv_y = info.emv.x(:,2,itraj,icell);
t = 0 : info.eyemvnt.conversionFactor : (length(emv_x)-1)*info.eyemvnt.conversionFactor;
img_size = info.img.imsize;

femc_plot_grat_sf_firing(grat_img,fr_spat,fr_temp,emv_x,emv_y,sf_list,info,t,nsf,img_size,ymax)

%% plot 2 - sensitivity v1 - by area under the curve
sscale = 100;
femc_plot_grat_sf_tuning(auc_spat,auc_temp,sf_list,info,sscale)

%% P cells simulation
fr_spat = zeros(ntpoint,2,nsf);
fr_temp = zeros(ntpoint,2,nsf); 
auc_spat = zeros(2,nsf);
auc_temp = zeros(2,nsf);
max_spat = zeros(2,nsf);
max_temp = zeros(2,nsf);

info.cell.type = 'P'; % P or M cell
for itest = 1:nsf
    
    info.img.stimsf = sf_list(itest);
    
    % Simulate LGN activity during letter viewing (2019.08.19)
    % This simulation contains 5 parts
    %  1. Stimulus by femc_stimulus.m
    %  2. Spatial Filtering
    %  3. Eye movement (imcomplete: need to work on 'simulated_BW' option (2018.09.25))
    %  4. Ramp protocol for contrast
    %  5. Temporal Filtering
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
    grat_img(:,:,:,itest) = img; 

    % TS = zeros(info.img.nimg,info.eyemvnt.nsamp,info.img.nstim);
    % pcproj = zeros(info.img.nimg,info.disc.pc,info.img.nstim,info.eyemvnt.ntraj,info.cell.ncell);
    % for icell = 1:info.cell.ncell
    %     for itraj = 1:info.eyemvnt.ntraj
    icell = 1;
    itraj = 1;
    for iimg = 1:info.img.nimg
        [ts_spatial] = femc_spatial(info,fimg(:,:,:,:,iimg),info.emv.x(:,:,itraj,icell)); 
%       [ts_ramp,ramp_contrast] = femc_ramp(ts_spatial,info);
        [Timeseries,info] = femc_temporal(model,ts_spatial,info);
%       [Timeseries,info] = femc_temporal(model,ts_ramp,info);
        TS(iimg,:,:) = Timeseries;
    end

    % plot the eye trajectories on the letters & firing rates
    img_size = info.img.imsize;
    emv_x = info.emv.x(:,1,itraj,icell);
    emv_y = info.emv.x(:,2,itraj,icell);
    t = 0 : info.eyemvnt.conversionFactor : (length(emv_x)-1)*info.eyemvnt.conversionFactor;
    ts_sfilt = squeeze(ts_spatial(:,1,:) - ts_spatial(:,2,:)); % center - surround
    ts_sfilt(ts_sfilt<0) = 0;
    % ts_max = round(max(Timeseries(:)),2);

    % save area under the curve 
    fr_spat(:,:,itest) = ts_sfilt;
    fr_temp(:,:,itest) = Timeseries;
    auc_spat(:,itest) = trapz(ts_sfilt);
    auc_temp(:,itest) = trapz(Timeseries);
    max_spat(:,itest) = max(ts_sfilt);
    max_temp(:,itest) = max(Timeseries);
    disp([num2str(sf_list(itest)) 'cpd done'])
end

%% plot 1 - sf vs firing 
ymax = [1.5*10^-4,0.01];
t = 0 : info.eyemvnt.conversionFactor : (length(emv_x)-1)*info.eyemvnt.conversionFactor;
img_size = info.img.imsize;
emv_x = info.emv.x(:,1,itraj,icell);
emv_y = info.emv.x(:,2,itraj,icell);

femc_plot_grat_sf_firing(grat_img,fr_spat,fr_temp,emv_x,emv_y,sf_list,info,t,nsf,img_size,ymax)

%% plot 2 - sensitivity v1 - by area under the curve
sscale = 25;
femc_plot_grat_sf_tuning(auc_spat,auc_temp,sf_list,info,sscale)
