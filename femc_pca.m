function [pcproj] = femc_pca(TS,info)
% reduce dimension of timeseries to principal components
            
% final version
[nimg,nt,nstim] = size(TS); % number of images (= n.samp = n.noise)
data = [];
for istim = 1:nstim
    data = [data;TS(:,:,istim)];
end

% PCA data in (row,column) = (n_observations, n_variables (or n_features))
pcproj = zeros(nimg,info.disc.pc,nstim);
coeff = pca(data);
% [coeff,score,latent,tsquared,explained,mu] = pca(data);
pcproj = [];
for istim = 1:nstim
    pcproj(:,:,istim) = squeeze(TS(:,:,istim)) * coeff(:,1:info.disc.pc);
end

% On 4/13/2022 - changed again 
% [U,S,V] = svd(data);
% reshape_TS = reshape(data,[nimg,nstim,nt]);
% for istim = 1:nstim
%     Vobs = V(:,1:info.disc.pc)' * squeeze(reshape_TS(:,istim,:))';
%     pcproj(:,:,istim) = Vobs';
%     
%     % Changed on May 5, 2020. I belive the old way (below) was wrong
%     % pcproj(:,:,istim) = squeeze(reshape_TS(:,istim,:)) * V(1:info.disc.pc,:)';
% end

% 2020/04/13 to be deleted. gettting rid of recenter
% recenter_data = data - repmat(mean(data),nimg*nstim,1);
% [~,~,V] = svd(recenter_data);
% reshape_rctTS = reshape(recenter_data,[nimg,nstim,nt]);
% for istim = 1:nstim
%     pcproj(:,:,istim) = squeeze(reshape_rctTS(:,istim,:)) * V(1:info.disc.pc,:)';
% end

% version 1
% n_train = n.img * n.train_percent;
% n_test = n.img - n_train;
% pcproj.train = zeros(n_train,n.pc,n.stim,n.traj,n.cell);
% pcproj.test  = zeros(n_test,n.pc,n.stim,n.traj,n.cell);
% kfold = round(1/(1-n.train_percent)); % k-fold validation (claculate from n.train_percent)
% rsTS = reshape(TS,[n.img/kfold,kfold,n.t,n.stim,n.traj,n.cell]);
% pcproj = zeros(n.img,n.pc,n.stim);
% for icell = 1:n.cell
%    for itraj = 1:n.traj
%         data = [];
%         for istim = 1:n.stim
%             data = [data;TS(:,:,istim)];
%         end
%         recenter_data = data - repmat(mean(data),n.img*n.stim,1);
%         [~,~,V] = svd(recenter_data);
%         reshape_rctTS = reshape(recenter_data,[n.img,n.stim,n.t]);
%         for istim = 1:n.stim
%             pcproj(:,:,istim) = squeeze(reshape_rctTS(:,istim,:)) * V(1:n.pc,:)';
%         end   
%    end
% end

% verson 2.
%         % seperate into training and testing sets for both stimuli
%         % row - observations; column - varaiables
%         dtrain = reshape(rsTS(:,1:kfold-1,:,:,itraj,icell),[n.img/kfold*(kfold-1),n.t,n.stim]);
%         dtest =  reshape(rsTS(:,kfold,:,:,itraj,icell),[n.img/kfold,n.t,n.stim]);       
%         d_train=[];
%         d_test=[];
%         for istim = 1:n.stim
%             d_train = [d_train;dtrain(:,:,istim)];
%             d_test = [d_test;dtest(:,:,istim)];
%         end
%         
%         % recenter the dataset
%         recenter_d_train = d_train - repmat(mean(d_train),n_train*n.stim,1);
%         recenter_d_test = d_test - repmat(mean(d_test),n_test*n.stim,1);
%         
%         % find the principal components of the training data
%         % coeff_train = pca(d_train);
%         % coeff_test = pca(d_test);
%         [~,~,V1] = svd(recenter_d_train);
%         [~,~,V2] = svd(recenter_d_test);
% 
%         % projections onto PCs
%         for istim = 1:n.stim
%             pcproj.train(:,:,istim,itraj,icell) = recenter_d_train(1+((istim-1)*n_train):istim*n_train,:) * V1(1:n.pc,:)';
%             pcproj.test(:,:,istim,itraj,icell)  = recenter_d_test(1+((istim-1)*n_test):istim*n_test,:) * V2(1:n.pc,:)';
%         end

end