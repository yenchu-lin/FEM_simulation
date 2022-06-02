function [img,noiseMask] = femc_Gaussian_noise(I,n)
% Adding Gaussian noise to the bitmap

% 1. create the Gaussian noise mask (default mu = 0)
mu = 0;
noiseMask = mu + sqrt(n.varnoise) * randn(size(I));

% 2. set the cutoff boundaries for Gaussian noise (boundary = cutoff*sd)
cutoff_boundary_hi = mu + n.cutoffnoise * sqrt(n.varnoise); 
cutoff_boundary_lo = mu - n.cutoffnoise * sqrt(n.varnoise); 
noiseMask = min(noiseMask, cutoff_boundary_hi);
noiseMask = max(noiseMask, cutoff_boundary_lo);

% 3. add the noise to the image
img = I + n.cnoise * noiseMask;

% 4. trim the signal if the noise mask range is outside of -1 to +1
img = min(img, 1);
img = max(img, -1);

% img = imcomplement(img);

end