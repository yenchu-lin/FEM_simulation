function [fd] = femc_fisherdisc(train_data,info)
% Use Fisher discriminant software to find the discrminant and the projections
% kfold = round(1/(1-n.train_percent)); % k-fold validation (claculate from n.train_percent)
% cd ../fisher_software
ntrain = info.img.nimg * info.disc.train_percent; 
ncell = info.cell.ncell;
nstim = info.img.nstim;
ntraj = info.eyemvnt.ntraj;

opts = fisherdisc_def();
opts.xvalid = 0; 
% 1 to do drop-one cross-validation and jackknife stats

opts.sub1 = 0; % default 0
% 1 to subtract 1 for empirical estimate of within-class covariances
% Note that this only matters if the number of elements in each class is unequal

opts.classifiers = {'mapequal'}; 
% opts.classifiers = {'mapbayes'}; 
% 'halfway' (halfway between the two midpoints)
% 'mapequal' (maximum a posteriori based on equal class occupancies)
% 'mapbayes' (maximum a posteriori based on empiric class occupancies)
% fd = cell(ncell,nstim);
for icell = 1:ncell % the cell number in one LGN simulation
    for itraj = 1:ntraj
        fddata = [];
                    
        % Rearrange the data into Fisher discrimination
        for istim = 1:nstim
            % row - feature (varaiables); 
            % column - sample (observations)
            fddata = [fddata,train_data(:,:,istim,itraj,icell)'];
        end
        
        % Iterate through the number of stimuli to run Fisher discriminant
        for ifd = 1:nstim
            tags = ones(1,ntrain*nstim)*2;
            tags(1,(ifd-1)*ntrain+1:ifd*ntrain) = 1;    
            [fd{icell,ifd}.results(itraj),fd{icell,ifd}.opts_used(itraj)] = fisherdisc(fddata,tags,opts);
%             [fd2{ifd}.results,fd2{ifd}.opts_used] = fisherdisc(fddata,tags,opts);
%             fd{icell,ifd}.results(itraj).discriminant = fd2{ifd}.results.discriminant;
        end
        
%         % Estimate the mean and variance of the a posteriori response ditributions
%         for ifd = 1:nstim
%             for istim = 1:nstim
%                 mu = fd2{ifd}.results.discriminant * fd2{ifd}.results.class_mean;
%                 % binary classifier -- class 1: target; class 2: non-target
%                 % mu and var of each distribution
%                 if ifd == istim
%                     iclass = 1;
%                 else 
%                     iclass = 2;
%                 end
%                 aprd.mu (iclass,ifd,icell,itraj) = mu(iclass);
%                 aprd.var(iclass,ifd,icell,itraj) = fd2{ifd}.results.var(iclass);
%             end
%         end
    end
end
end