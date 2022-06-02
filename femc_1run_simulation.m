%% Simulate LGN activity and behavior responses during letter viewing
% Last update on 2022/04/24 by YCL 

%% 0. make letters (jpg to mat)
% import jpg from illustrator
% by femc_make_bmp_letters.m
% make a mat file for each letter (512 x 512)
% store in folder - letter_mat

%% Inputs
newstructname = 'Test2_5M_Gaussian_JI_c1.mat';
savedir = '/Users/mb/Documents/Data/Femc_simulation/results_092121/';
disp(newstructname)

% 1. input stimuli 
% 1.1. input type - letters
info.img.type = 'letter';
info.img.nstim = 2;
info.stim = cell(info.img.nstim,1);
info.stim{1} = 'H';
info.stim{2} = 'N';
% info.stim{1} = 'Eup';
% info.stim{2} = 'Edown';
% info.stim{3} = 'Eright';
% info.stim{4} = 'Eleft';
% --------------------------
% 1.2. input type - gratings 
% info.img.type = 'grating';
% info.img.nstim = 2;
% info.img.stimsf = 4;
% info.stim = cell(info.img.nstim,1);
% info.stim{1} = 'vertical';
% info.stim{2} = 'horizontal';
% info.img.shift = 0; % shift the amount of the spatial period of the stimulus; one-quarter (0.25) or one-half (0.5)

% 2.image parameters
info.img.imsize = 512;
info.img.D = 2;
info.img.npad = 2;
info.img.nimg = 100;

% 3. noise parameters
info.img.noisetype = 'Gaussian'; % 1/F or Gaussian
info.img.cnoise = 0.3;
info.img.cstim  = 0.01; % 0.0003
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
info.eyemvnt.type = 'exp_traj'; % 'exp_traj' or 'sim_traj'
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
info.cell.type = 'M'; % P or M cell
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

% Interactive Inputs
% info.img.nstim = getinp('how many stimuli (1-26 letters)','d',[1,26],4);
% for istim = 1:info.img.nstim
%     info.stim{istim} = getinp('input letter (A to Z)','s',[],'E');
% end
% info.img.imsize = getinp('number of pixels per image array','d',[64,2048],128);
% info.img.D = getinp('stimulus display size (in deg)','d',[1,10],2); 
% info.img.cstim = getinp('contrast of the stimulus (-1 to 1)','f',[-1,1],0.04);
% info.img.npad = 2; % add zero padding based on npad (default npad = 2)
% info.img.nimg = getinp('number of trials','d',[1,1000],100);
% info.img.cnoise = getinp('contrast of the noise  (-1 to 1)','f',[-1,1],0.3);
% info.img.varnoise = getinp('variance of the Gaussian noise','f',[0,1],1);
% info.img.cutoffnoise = getinp('how many folds of std as the cutoff range for Gaussian noise','d',[1,4],2);
% info.eyemvnt.type = getinp('input eye movement','s',[],'exp_traj'); % trajactory from single trial (experimental data)
% info.cell.temp.delta_t = getinp('the time step (sec)?','f',[0,1],0.001); % sec; the time step for filtering 

% Simulate LGN activity during letter viewing (2019.08.19)
% This simulation contains 5 parts
%  1. Stimulus by femc_stimulus.m
%  2. Spatial Filtering
%  3. Eye movement (imcomplete: need to work on 'simulated_BW' option (2018.09.25))
%  4. Ramp protocol for contrast
%  5. Temporal Filtering
[model,info] = femc_lgn_parameters(info);
[info] = femc_cell_center(info,model);
% [info] = femc_cell_center_old(info,model); % for testing purpose
[info] = femc_input_eye_movement(info);

% istim/itraj/icell/iimg
imsize = info.img.imsize * info.img.npad;
fimg = zeros(imsize,imsize,info.img.nstim,2,info.img.nimg);
for iimg = 1:info.img.nimg
    [img,img_rotate,noiseMask] = femc_stimulus(info);
    [filterimg] = femc_spatial_filter_img(info,img,model);
    fimg(:,:,:,:,iimg) = filterimg;
end
disp('filter img done')

info.eyemvnt.ntraj = 10;
TS = zeros(info.img.nimg,info.eyemvnt.nsamp,info.img.nstim);
pcproj = zeros(info.img.nimg,info.disc.pc,info.img.nstim,info.eyemvnt.ntraj,info.cell.ncell);
for icell = 1:info.cell.ncell
    for itraj = 1:info.eyemvnt.ntraj
        for iimg = 1:info.img.nimg
            [ts_spatial] = femc_spatial(info,fimg(:,:,:,:,iimg),info.emv.x(:,:,itraj,icell)); 
            [ts_ramp,ramp_contrast] = femc_ramp(ts_spatial,info);
            [Timeseries,info] = femc_temporal(model,ts_ramp,info);
            TS(iimg,:,:) = Timeseries;
        end
        
        % 4. PCA
        % pc_value - results of principal component projections
        % pcproj(n.nimg, n.pc, n.stim, n.traj, n.cell)
        [proj] = femc_pca(TS,info);
        pcproj(:,:,:,itraj,icell) = proj;
    end
    disp('cell')
    disp(icell)
end

% Decision model
% This simulation contains 3 parts
% N-Fold Validation
% pcproj(n.nimg, n.pc, n.stim, n.traj, n.cell)
% n.train_percent = 0.80;
nfold = round(1/(1-info.disc.train_percent));
ntest = info.img.nimg * (1-info.disc.train_percent);
ntrain = info.img.nimg * info.disc.train_percent;
fc_test_pop  = zeros(info.eyemvnt.ntraj, info.cell.ncell, nfold);
fc_test_ind  = zeros(info.eyemvnt.ntraj, info.cell.ncell, nfold);
for ifold = 1:nfold
    
    % 5. Fisher discriminator by training data
    % fd{icell,ifd}.results(itraj) 
    % separate training and testing data
    ind1 = round(1+(ifold-1)*ntest);
    ind2 = round(ifold*ntest);
    test_data = pcproj(ind1:ind2,:,:,:,:);
    train_data = pcproj;
    train_data(ind1:ind2,:,:,:,:)=[];
    
    [fd] = femc_fisherdisc(train_data,info);
    disp([num2str(ifold), ' fd done'])

    % Estimate the mean and variance of the Gaussians from each cell
    % aprd(stim*2, cell, traj)
    cd ../eye_movement_control_simulation
    [aprd] = femc_aprd(fd,info);
    disp([num2str(ifold), 'aprd done'])

    % 6. Calculate the confusion matrix and fraction correct
    fc = zeros(info.eyemvnt.ntraj,info.cell.ncell);
    cell_no = [];
    for icell = 1:info.cell.ncell
    
        % Calculate the log likelihood by training set & testing set
        % Estimate the log likelihood ratio of all cells
        cell_no = [cell_no, icell];
        [cm_test_pop] = femc_loglikelihood(test_data,fd,aprd,info,cell_no);
        [cm_test_ind] = femc_loglikelihood(test_data,fd,aprd,info,icell);

        % Estimate the fraction correct based on cm
        [fcpop] = femc_fraction_correct(cm_test_pop,info);
        [fcind] = femc_fraction_correct(cm_test_ind,info);
        fc_test_pop(:,icell,ifold) = fcpop;
        fc_test_ind(:,icell,ifold) = fcind;
    end
    disp([num2str(ifold), 'fc done'])
end

% Save the result
data.info = info;
data.img = img;
data.img_rotate = img_rotate;
data.pcproj = pcproj;
data.fd = fd;
data.aprd = aprd;
data.fc.pop = fc_test_pop;
data.fc.ind = fc_test_ind;

% save the results
savename = newstructname;
save(fullfile(savedir,savename),'data', '-v7.3')
% save(fullfile(savedir,savename),'data', 'Variablename', '-v7.3')
% save(newstructname,'-struct','data')
fprintf(1, '\n');
disp([savename ' file saved'])
