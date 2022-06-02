function [ts_spatial] = femc_spatial(info,filterimg,emv)
% The function move the eye trajectory on the spatial filtering images
    % currently there is only one type option - circular Gaussian
% Also see: femc_lgn_spatial_circular_Gaussian.m
%           femc_spatial_filter.m
% Last updated on 2022/04/04

imsize = info.img.imsize * info.img.npad;
ts_spatial = zeros(length(emv),2,info.img.nstim); % 2 is for iteration of receptive fields: center and surround 

% Round the coordinates (can be improved) 
% roundx = round(emv.x * (imsize/2)); % for simulated eye traj
roundx = round(emv);

% To deal with the eye movemement outside of the image, set
% the coordinates as (imsize,imsize)
roundx((roundx>imsize)) = imsize; 
roundx((roundx<=0)) = imsize;
            
for irf = 1:2
    for istim = 1:info.img.nstim
        
%         Img = flipud(filterimg(:,:,istim,irf)); % updated on 2022/04/04  
        Img = filterimg(:,:,istim,irf);  
        img_size = size(Img);
        index = (roundx(:,1)-1)*img_size(1) + roundx(:,2);
        allpixels = reshape(Img,1,img_size(1)*img_size(2));
        ts_spatial(:,irf,istim) = allpixels(index)';

    end
end

end