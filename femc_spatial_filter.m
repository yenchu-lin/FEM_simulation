function [filterImg] = femc_spatial_filter(img,g)
% recreate spatial filtered image by FT back

ftImg = fft2(img);
ftFilterImg = ftImg .* g; % multiply under freq domain to avoid convolution
filterImg = ifft2(ftFilterImg);

end