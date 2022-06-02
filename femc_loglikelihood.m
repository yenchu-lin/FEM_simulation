function [cm_test] = femc_loglikelihood(test_data,fd,aprd,info,tbl)
% Estimate the log likelihood of each trajectory
% by adding the likelihood across cells

nimg = info.img.nimg;
nstim = info.img.nstim;
ncell = info.cell.ncell;
ntraj = info.eyemvnt.ntraj;

% Responses from all stimuli
ntest = round(nimg * (1-info.disc.train_percent));
proj_test = zeros(ntest,nstim,nstim,ncell,ntraj);
for itraj = 1:ntraj
    for icell = 1:ncell
        for istim = 1:nstim
            for ifd = 1:nstim
                proj_test(:,istim,ifd,icell,itraj)  = test_data(:,:,istim,itraj,icell) * fd{icell,ifd}.results(itraj).discriminant';
            end
        end
    end
end

% Estimate the log likelihood of each response
% log P(r|target) ; stim  = fd ---> class 1
% log P(r|~target); stim ~= fd ---> class 2
ll_test  = zeros(ntest,nstim,nstim,ntraj,ncell,2);
for icell = 1:ncell 
    for itraj = 1:ntraj
        for ifd = 1:nstim
            for istim = 1:nstim
                for iclass = 1:2 % binary classification
                    mu = aprd.mu(iclass,ifd,icell,itraj);
                    v = aprd.var(iclass,ifd,icell,itraj); 
                    ll_test(:,istim,ifd,itraj,icell,iclass) = -((proj_test(:,istim,ifd,icell,itraj)-mu).^2 / (2*v)) - 0.5*log(2*pi*v);
                end
            end
        end
    end
end

% Estimate the log likelihood ratio of each cell
% output: log_likelihood_ratio(response, traj, cell)
llr_test  = ll_test(:,:,:,:,:,1) -  ll_test(:,:,:,:,:,2); % substract the 6th dimension - P(target) True vs. P(non-target)

% Combine all the cells from the selected subset based on the subset table
% output: combined (total) log_likelihood_ratio(response, stim, traj)
tllr_test  = sum(llr_test(:,:,:,:,tbl),5); % add across the 5th dimension - cells

% Confusion matrix
cm_test  = zeros(nstim,nstim,ntraj);
for itraj = 1:ntraj
    for irow = 1:nstim % ifd
        [~,ind_test] =  max(tllr_test(:,:,irow,itraj),[],2);
        for icol = 1:nstim % istim
            cm_test(irow,icol,itraj)  = sum(ind_test == icol);
        end
    end
end

end