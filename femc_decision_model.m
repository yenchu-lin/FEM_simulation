function [fc_test_pop,fc_test_ind] = femc_decision_model(pcproj,info)
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
%     [fd] = femc_fisherdisc(train_data,info);
%     disp([num2str(ifold), ' fd done'])
% 
    % Estimate the mean and variance of the Gaussians from each cell
    % aprd(stim*2, cell, traj)
    [aprd] = femc_aprd(fd,info);
%     disp([num2str(ifold), 'aprd done'])

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
end
