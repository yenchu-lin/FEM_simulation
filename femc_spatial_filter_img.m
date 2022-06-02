function [filterimg] = femc_spatial_filter_img(info, img, model)
% The function process spatial filtering
    % currently there is only one type option - circular Gaussian
% Also see: femc_lgn_spatial_circular_Gaussian.m

imsize = info.img.imsize * info.img.npad;
filterimg = zeros(imsize, imsize, info.img.nstim, 2);
for istim = 1:info.img.nstim
        
        % Center; RF in freq domain
        RF = 'center';
        [~,g] = femc_lgn_spatial_circular_Gaussian(RF,info,model); % g is in freq domain; range from 0 to 1 
        [filterImg_C] = femc_spatial_filter(img(:,:,istim),g);
        filterimg(:,:,istim,1) = filterImg_C;

        % Surround; RF in freq domain
        RF = 'surround';
        [~,g] = femc_lgn_spatial_circular_Gaussian(RF,info,model); % g is in freq domain; range from 0 to 1 
        [filterImg_S] = femc_spatial_filter(img(:,:,istim),g);
        filterimg(:,:,istim,2) = filterImg_S;
        
end

end