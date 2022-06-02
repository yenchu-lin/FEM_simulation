% Test variuos numbers of PCs vs. performance
% Simulate LGN activity and behavior responses during letter viewing
% Last update on 2022/04/13 by YCL 

% Inputs
%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%
disp('Test the weights of individual PCs')
npc = 10; % number of pc tested

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
info.img.cnoise = 0.3; % 0.3;
info.img.cstim  = 0.0001; % 0.0003
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
info.cell.type = 'M'; % P or M cell
info.cell.rftype = 'ON'; % ON or OFF cell
info.cell.spat.model_no = 1; % 1- circular Gaussian
info.cell.temp.delta_t = 0.001;
info.cell.temp.model_no = 3; % 2- gamma; 3- HPLP

% 6. discriminator parameters
info.disc.train_percent = 0.80; % percentage of data into the training 
info.disc.pc = 10; % number of principal components
%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%

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

% test using different number of PCs
info.cell.ncell = 5;
info.eyemvnt.ntraj = 10;
fc_test_pop = cell(npc,1);
fc_test_ind = cell(npc,1);
for ipc = 1:npc
    pc_proj= zeros(info.img.nimg, ipc, info.img.nstim, info.eyemvnt.ntraj, info.cell.ncell);
    TS = zeros(info.img.nimg,info.eyemvnt.nsamp,info.img.nstim);
    for icell = 1:info.cell.ncell
        for itraj = 1:info.eyemvnt.ntraj
            for iimg = 1:info.img.nimg
                [ts_spatial] = femc_spatial(info,fimg(:,:,:,:,iimg),info.emv.x(:,:,itraj,icell)); 
                [ts_ramp,ramp_contrast] = femc_ramp(ts_spatial,info);
                [Timeseries,info] = femc_temporal(model,ts_ramp,info);
                TS(iimg,:,:) = Timeseries;
            end

            % PCA
            [nimg,nt,nstim] = size(TS);
            data = []; 
            for istim = 1:nstim
                data = [data;TS(:,:,istim)]; 
            end

            % PCA data in (row,column) = (n_observations, n_variables (or n_features))
            [coeff,score,latent,tsquared,explained,mu] = pca(data);
            pcproj = [];
            for istim = 1:nstim
                pcproj(:,:,istim) = squeeze(TS(:,:,istim)) * coeff(:,1:ipc);
            end
            pc_proj(:,:,:,itraj,icell) = pcproj;
        end            
    end

    % decision model - fraction correct of single cell and single trajactory
    info.disc.pc = ipc;
    [fc_test_pop{ipc},fc_test_ind{ipc}] = femc_decision_model(pc_proj,info);
end

%% plot - look at performance vs. number of PCs
% plot pc performance of a single cell
% fc_test_ind{ipc}(itraj,icell,ifold)
saveplot = 0;
itraj = 1;
icell = 1;
% average across the 5 fold validation
fc_ind = [];
for ipc = 1:npc
    fc_ind(ipc,1) = mean(fc_test_ind{ipc}(itraj,icell,:));
end

% plot 1 - number of pc vs. performance
figure
plot(1:npc,fc_ind,'o-','linewidth',2)
xlabel('number of pc used')
ylabel('performance')
set(gca,'fontsize',18,'linewidth',2)
box off
title([info.cell.type ' ' info.cell.rftype ' cell (1 cell, 1 traj)'])

% save figure
fname = ['ch2_pc_performance' info.cell.type ' cell'];
figpt = '/Users/mb/Desktop/Thesis/Figures';
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end

%% plot - look at performance vs. number of PCs
% plot pc performance of multiple cells
saveplot = 0;
% average across the 5 fold validation
fc_ind = [];
for ipc = 1:npc
    for icell = 1:info.cell.ncell
        fc_pop(icell,ipc) = mean(fc_test_pop{ipc}(itraj,icell,:));
    end
end

% plot - number of pc vs. performance
figure('Position',[10 10,800,600])
hold on
for icell = 1:info.cell.ncell
    plot(1:npc,fc_pop(icell,:),'o-','linewidth',2)
end
xlabel('number of pc used')
ylabel('performance')
ylim([0.4 1])
legend('1 cell','2 cells','3 cells','4 cells','5 cells','Location','southeast')
set(gca,'fontsize',18,'linewidth',2)
box off
title([info.cell.type ' ' info.cell.rftype ' cell'])

% save figure
fname = ['ch2_pc_performance' info.cell.type ' cell_multicells'];
figpt = '/Users/mb/Desktop/Thesis/Figures';
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end

%% plot - scatter plot of two pcs
pc1 = 1;
pc2 = 2;
iimg = 1;
figure
hold on
istim = 1; % info.stim{istim}
plot(pcproj(:,pc1,istim), pcproj(:,pc2,istim),'bo')
istim = 2; % info.stim{istim}
plot(pcproj(:,pc1,istim), pcproj(:,pc2,istim),'ro')
% plot(0:0.01:0.15,0:0.01:0.15,'--','color',[.7 .7 .7])
legend(info.stim{1},info.stim{2})
xlabel(['pc' num2str(pc1)])
ylabel(['pc' num2str(pc2)])
set(gca,'fontsize',18,'linewidth',2)
box off
title([info.cell.type ' ' info.cell.rftype ' cell (1 cell, 1 traj)'])
axis equal 
axis square

%% sanity check for PCA 
% img_size = info.img.imsize;
% emv_x = info.emv.x(:,1,itraj,icell);
% emv_y = info.emv.x(:,2,itraj,icell);
% t = 0 : info.eyemvnt.conversionFactor : (length(emv_x)-1)*info.eyemvnt.conversionFactor;
% ts_sfilt = squeeze(ts_spatial(:,1,:) - ts_spatial(:,2,:)); % center - surround

% [coeff,score,latent,tsquared,explained,mu] = pca(data);
% [U,S,V] = svd(data);

% Calculate eigenvalues (eigV) and eigenvectors (eigD) of the covariance matrix
% covarianceMatrix = cov(data);
% [eigV,eigD] = eig(covarianceMatrix);
% [coeff_svd,latent_svd] = svd(covarianceMatrix); % Singular value decomposition

% Multiply the original data by the principal component vectors
% - to get the projections of the original data on the principal component 
% vector space. This is also the output "score"
% reshape_TS = reshape(data,[nimg,nstim,nt]);
% for istim = 1:nstim
%     Vobs = V(:,1:info.disc.pc)' * squeeze(reshape_TS(:,istim,:))';
%     Vobs = V' * squeeze(reshape_TS(:,istim,:))';
%     pcproj(:,:,istim) = Vobs';
%     pcproj(:,:,istim) = squeeze(TS(:,:,istim)) * coeff;
% end
