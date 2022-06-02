function [img,noiseMask] = femc_OneOverF_noise2(I,imginfo)
% Add 1/f noise to the bitmap

% Create random noise
imsize = round(size(I,2));
im_rand = rand(imsize);
freup = imginfo.noisefrequp;
frelo = imginfo.noisefreqlo;

% Builds X and Y of distance from centre
stepsize = (freup-frelo)/(imsize/2);
[X,Y] = meshgrid([ freup-stepsize:-stepsize:frelo, frelo:stepsize:freup-stepsize ],...
                 [ freup-stepsize:-stepsize:frelo, frelo:stepsize:freup-stepsize ]);
% stepsize = (freq_limit_up-1)/(imsize/2);             
% [X,Y] = meshgrid([ freq_limit_up-stepsize:-stepsize:1, 1:stepsize:freq_limit_up-stepsize ],...
%                  [ freq_limit_up-stepsize:-stepsize:1, 1:stepsize:freq_limit_up-stepsize ]);
% [X,Y] = meshgrid([ imsize/2:-1:1, 1:imsize/2 ],...
%                  [ imsize/2:-1:1, 1:imsize/2 ]);

% Creates radial matrix
R = sqrt(X.^2 + Y.^2);

% FFT of image and matrix inversion
IM1 = fftshift(fft2(im_rand));
IM2 = IM1./abs(IM1); % Set magnitude to one

% Apply 1/F.^2
F = 1./R;
IM3 = IM2 .* F;  

% Transform to space domain
n = real( ifft2(ifftshift(IM3)) );
nm = n ./ sqrt(sum(n(:).^2)) ; % normalize by the square root of the sum of squares
noiseMask = nm;

% add the noise to the image
img = I + imginfo.cnoise * noiseMask;
% img = I + noiseMask;

% trim the signal if the noise mask range is outside of -1 to +1
img = min(img, 1);
img = max(img, -1);

% [m,n] = size(I);
% f = imginfo.freq;      % frequency (in Hz) of 1/f^2 noise mask 
% m = round(m); 
% n = round(n); 
% N = m*n;
% x = randn(1, N); % generate white noise sequence
% X = fft(x);
% 
% % prepare a vector with frequency indexes 
% NumUniquePts = N/2 + 1;              % number of the unique fft points
% k = linspace(1,f,NumUniquePts);      % vector with frequency indexes 
% 
% % manipulate the left half of the spectrum so the PSD
% % is proportional to the frequency by a factor of 1/f, 
% % i.e. the amplitudes are proportional to 1/sqrt(f)
% X = X(1:NumUniquePts);      
% X = X./sqrt(k);
% 
% % prepare the right half of the spectrum - a conjugate copy of the left
% % one except the DC component and the Nyquist component - they are unique,
% % and reconstruct the whole spectrum
% Xlist = [X(1:end-1) real(X(end)) fliplr(conj(X(2:end-1)))];
% y = real(ifft(Xlist)); % IFFT
% 
% % ensure that the length of y is N
% y = y(1, 1:N);
% 
% % form the noise matrix and ensure unity standard 
% % deviation and zero mean value (columnwise)
% y = reshape(y, [m, n]);
% y = bsxfun(@minus, y, mean(y));
% y = bsxfun(@rdivide, y, std(y));

% figure
% imagesc(y)
% colormap(gray)
% figure
% loglog(k, abs(X)) % plot the magnitude of K
% xlabel('Frequency (Hz)')
% ylabel('Response (impulses/sec)')

end

