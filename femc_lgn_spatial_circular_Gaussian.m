function [F,scaled_g] = femc_lgn_spatial_circular_Gaussian(RF,info,model)
% spatial response profile by circular Gaussian
%   see reference Rucci, Edelman & Jonathan Wray, 2000, J. Neuroscience

% parameters for the LGN model
filtertype = 'spat';
[p] = femc_lgn_multiplier(RF,model,info,filtertype);

% p.strength -- amplitude
% p.radius -- respective size in degree

N = info.img.imsize * info.img.npad;

X=[];
X(1,:)=linspace(-info.img.D,info.img.D,N); % in deg
X(2,:)=linspace(-info.img.D,info.img.D,N); % in deg

% Gaussian RF in space domain
F = p.str/ (2*pi*(p.radius^2)) * exp((-(X(1,:).^2 + X(2,:).^2)) ...
    / (2*(p.radius ^ 2))); 

% create an array of omega (w)
wx = min((0:N-1), N-(0:N-1)); 
wy = wx';
warray = sqrt(repmat(wx.^2,N,1)+repmat(wy.^2,1,N)) * (2*pi/info.img.D/info.img.npad); % in radian/deg

% Gaussian RF in freq domain
g = p.str * exp(-(warray.^2) * (p.radius^2) / 2); % the values are all positive
scaled_g = g / max(g(:)); % rescale to the range from 0 to 1

% Plot
% ticksize = D*2; % number of x and y axis labeling ticks (in deg) for plotting
% figure
% imagesc(warray)
% colorbar
% title('w array for the temporal profile of RF')

% figure
% subplot(1,2,1)
% x = linspace(0,D/2,N/2); % mod(n,2) == 0 two zeros in the middle (need fix!)
% x = [fliplr(-x),x];
% plot(x,F)
% title('spatial Gaussian RF')
% set(gca, 'FontSize', 16)
% subplot(1,2,2)
% imagesc(log(g))
% xlabel('x (deg)')
% 
% title('RF in freq domain (log scale)')
% set(gca, 'FontSize', 16)
% A=(-(D*npad/2):D*npad/ticksize:D*npad/2)';
% s = num2str(A);
% set(gca,'XTick',(0:1/ticksize:1)*N,'XTickLabel',s)
% set(gca,'YTick',(0:1/ticksize:1)*N,'YTickLabel',s)
% xlabel('x (deg)')
% ylabel('y (deg)')

% figure
% imagesc(log(scaled_g))
% xlabel('x (deg)')
% title('RF in freq domain (log scale)')
% set(gca, 'FontSize', 16)
% A=(-(D*npad/2):D*npad/ticksize:D*npad/2)';
% s = num2str(A);
% set(gca,'XTick',(0:1/ticksize:1)*N,'XTickLabel',s)
% set(gca,'YTick',(0:1/ticksize:1)*N,'YTickLabel',s)
% xlabel('x (deg)')
% ylabel('y (deg)')

end