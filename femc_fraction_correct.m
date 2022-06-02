function [accur] = femc_fraction_correct(cm,info)
% Estimate the fraction correct based on the confusion matrix
% Reference: https://www.dataschool.io/simple-guide-to-confusion-matrix-terminology/
ntraj = info.eyemvnt.ntraj;

% Accuracy
accur = zeros(ntraj,1);
for itraj = 1:ntraj
    accur(itraj,1) = sum(diag(cm(:,:,itraj)))/sum(sum(cm(:,:,itraj)));
end

end