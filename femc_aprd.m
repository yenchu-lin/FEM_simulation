function [aprd] = femc_aprd(fd,info)
% Estimate the mean and variance of the a posteriori response ditributions
% modified version (April 23, 2019)

ncell = info.cell.ncell;
ntraj = info.eyemvnt.ntraj;
nstim = info.img.nstim;

% reorganize fd data to get the projections of each condition from FD
% ntrain = round(n.trial * n.train_percent);
aprd.mu =  zeros(2,nstim,ncell,ntraj);
aprd.var = zeros(2,nstim,ncell,ntraj);
for icell = 1:ncell
    for itraj = 1:ntraj
        for ifd = 1:nstim
            for istim = 1:nstim
                mu = fd{icell,ifd}.results(itraj).discriminant * fd{icell,ifd}.results(itraj).class_mean;
                % binary classifier -- class 1: target; class 2: non-target
                % mu and var of each distribution
                if ifd == istim
                    iclass = 1;
                else 
                    iclass = 2;
                end
                aprd.mu (iclass,ifd,icell,itraj) = mu(iclass);
                aprd.var(iclass,ifd,icell,itraj) = fd{icell,ifd}.results(itraj).var(iclass);
            end
        end
    end
end

end