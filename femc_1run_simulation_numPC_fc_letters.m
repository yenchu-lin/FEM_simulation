% Test variuos numbers of PCs vs. performance while using two letter pairs
% HN and EF
% Note that all the 5 cells are at the center therefore perform worse for EF
%
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
info.img.cstim  = 0.0002; % for P cell use 0.0002; for M cell use 0.0001
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
    disp([num2str(ipc) ' pc done'])
end

% save the performance
fc.HN.fcpop = fc_test_pop;
fc.HN.fcind = fc_test_ind;
disp('HN done')
%
disp('letter EF')
%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%
info.stim = cell(info.img.nstim,1);
info.stim{1} = 'E';
info.stim{2} = 'F';
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
    
    disp([num2str(ipc) ' pc done'])
end

% save the performance
fc.EF.fcpop = fc_test_pop;
fc.EF.fcind = fc_test_ind;
disp('EF done')

% average across 5-fold validation
for ipc = 1:npc
    for itraj = 1:info.eyemvnt.ntraj
        for icell = 1:info.cell.ncell
            fc_pop.HN(icell,itraj,ipc) = mean(fc.HN.fcpop{ipc}(itraj,icell,:));
            fc_pop.EF(icell,itraj,ipc) = mean(fc.EF.fcpop{ipc}(itraj,icell,:));
        end
    end
end
% average across all 10 trajectories (performance of 5 cells together)
icell = 5;
for ipc = 1:npc
    for itraj = 1:info.eyemvnt.ntraj
        fc_mean.HN(ipc) = mean(fc_pop.HN(icell,:,ipc));
        fc_mean.EF(ipc) = mean(fc_pop.EF(icell,:,ipc));
    end
end

%% plot
saveplot = 1;
cmap = colormap(parula);
figure('Position',[10 10,800,600])
hold on
for ipc = 1:npc
    x = fc_mean.HN(ipc);
    y = fc_mean.EF(ipc);
    plot(x,y,'o','linewidth',2,'color',cmap(ipc*25,:))
    text(x+0.002,y+0.002,[num2str(ipc) ' pc'],'fontsize',16)
end

xlabel('performance (H vs. N)')
ylabel('performance (E vs. F)')
% xlim([.65 .75])
% ylim([.5 .6])
axis equal
axis square
legend('1 pc','2 pcs','3 pcs','4 pcs','5 pcs','6 pcs','7 pcs','8 pcs','9 pcs','10 pcs','Location','southeast')
set(gca,'fontsize',18,'linewidth',2)
box off
% title('dicrimination performance') % title([info.cell.type ' ' info.cell.rftype ' cell'])

% save figure
fname = ['ch2_pc_performance_' info.cell.type 'cell_letterHN_EF_sim2'];
figpt = '/Users/mb/Desktop/Thesis/Figures';
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end
