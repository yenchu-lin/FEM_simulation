%% Test P and M cell tuning curves
% Simulate spatial dn temporal filter sensitivity given differene freq (w)
% Last update on 2022/04/07 by YCL 

%% M or P cells simulation
%%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%
% neural population parameters
info.cell.ncell = 1; % number of cells
info.cell.type = 'M'; % P or M cell
info.cell.rftype = 'OFF'; % ON or OFF cell
info.cell.spat.model_no = 1; % 1- circular Gaussian
info.cell.temp.delta_t = 0.001;
info.cell.temp.model_no = 3; % 2- gamma; 3- HPLP
%%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%

[model,info] = femc_lgn_parameters(info);
[info] = femc_cell_center_old(info,model); 

% spatial filtering
% see paper 
% Receptive Fields of P and M Ganglion Cells Across the Primate Retina
% LISA J. CRONER, EHUD KAPLAN (1995)
%%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%
saveplot = 1;
filtertype = 'spat';
N = 1000;
X = linspace(-1,1,N);
%%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%

RF = 'center';
[pc] = femc_lgn_multiplier(RF,model,info,filtertype);
Rc = pc.Kc * pc.rc * exp(-(X.^2/pc.rc^2)); % response in space domain

% Fc = p.str * exp(- ((wspat-X).^2) / (p.radius^2) ); 

RF = 'surround';
[ps] = femc_lgn_multiplier(RF,model,info,filtertype);
Rs = ps.Ks * ps.rs * exp(-(X.^2/ps.rs^2)); % response in space domain

% these are the response in freq domain
% p.strength -- amplitude
% p.radius -- respective size in degree
% wspat = linspace(0,5,N); % freq list
% Rc = pc.Kc * pi * pc.rc^2 * exp(-(pi * pc.rc .* wspat)).^2;
% Rs = ps.Ks * pi * ps.rs^2 * exp(-(pi * ps.rs .* wspat)).^2;

% center - surround
R_sp = Rc - Rs;

% plots 
% plot spatial filter
% xplot = linspace(-5,5,N);
figure('Position',[2000 10 600 600])
hold on
plot(X,R_sp,'linewidth',2)
plot(X,Rc,'--','color',[.6 .6 .6],'linewidth',2)
plot(X,-Rs,'--','color',[.8 .8 .8],'linewidth',2)
legend('spatial filter','center','surround','Location','southeast')
xlabel('X (deg)')
title([info.cell.type, ' ', info.cell.rftype, ' cell'])
set(gca,'fontsize',18,'linewidth',2)
box off

% plot - save figure
fname = ['ch2_tuning_curve_' info.cell.type, '_', info.cell.rftype, '_cell'];
figpt = '/Users/mb/Desktop/Thesis/Figures';
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end

%% Temporal filtering
%%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%
saveplot = 0;
% N = 100;
% wtemp = logspace(0,2); % freq list
% wtemp = linspace(0,64,N); % freq list
wtick = [0.1,0.25,0.5,1,2,4,8,16,32,64,128];
N = 2^11;
nhalf = 1/2*N;
fmax = 64;
wmax = fmax*2*pi;
wmin = wmax/nhalf;
wtemp = (0:nhalf)*wmin;
%%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%

filtertype = 'temp';
RF = 'center';
[p] = femc_lgn_multiplier(RF,model,info,filtertype);
Kc = p.A * exp(-1i * wtemp * p.D) .* (1 - p.Hs ./ (1 + 1i * wtemp * p.tau_S))...
    .* (1 ./ (1 + 1i * wtemp * p.tau_L).^ p.NL); % radian/sec

RF = 'surround';
[p] = femc_lgn_multiplier(RF,model,info,filtertype);

% in freq domain
Ks = p.A * exp(-1i * wtemp * p.D) .* (1 - p.Hs ./ (1 + 1i * wtemp * p.tau_S))...
    .* (1 ./ (1 + 1i * wtemp * p.tau_L).^ p.NL); % radian/sec

wlist = abs(ifft(wtemp));

% plots 
% plot temporal filter
figure('Position',[2300 10 600 1000])
subplot(2,1,1)
% semilogx(wtemp/2/pi,abs(Kc),'linewidth',2)
loglog(wtemp/2/pi,abs(Kc),'linewidth',2)
xticks(wtick)
legend('center','Location','southeast')
title([info.cell.type, ' ', info.cell.rftype, ' cell'])
set(gca,'fontsize',18,'linewidth',2)
box off

subplot(2,1,2)
% semilogx(wtemp/2/pi,abs(Ks),'linewidth',2)
loglog(wtemp/2/pi,abs(Ks),'linewidth',2)
xticks(wtick)
xlabel('frequency (Hz)')
legend('surround','Location','southeast')
set(gca,'fontsize',18,'linewidth',2)
box off

% subplot(3,1,3)
% semilogx(wtemp/2/pi,abs(Kc-Ks),'linewidth',2)
% xticks(wtick)
% legend('temporal filter','Location','southeast')
% set(gca,'fontsize',18,'linewidth',2)
% box off

% plot - save figure
fname = ['ch2_tuning_curve_temporal_' info.cell.type, '_', info.cell.rftype, '_cell'];
figpt = '/Users/mb/Desktop/Thesis/Figures';
if saveplot
    saveas(gca,fullfile(figpt,fname), 'tif')
end

%% Spatial filtering
% see paper 
% "Spatiotemporal Receptive Field Organization in the Lateral Geniculate Nucleus of Cats and Kittens"
% (1997) DAQING CAI, GREGORY C. DEANGELIS, AND RALPH D. FREEMAN
% %%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%
% N = 1000;
% X = linspace(-1,4,N);
% wspat = linspace(0,3,N); % freq list
% saveplot = 0;
% %%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%
% 
% filtertype = 'spat';
% % p.strength -- amplitude
% % p.radius -- respective size in degree
% RF = 'center';
% [p] = femc_lgn_multiplier(RF,model,info,filtertype);
% Fc = p.str * exp(- ((wspat-X).^2) / (p.radius^2) ); 
% 
% RF = 'surround';
% [p] = femc_lgn_multiplier(RF,model,info,filtertype);
% Fs = p.str * exp(- ((wspat-X).^2) / (p.radius^2) ); 
% 
% F_sp = Fc - Fs;
% 
% % plots 
% % plot spatial filter
% xplot = linspace(0,3,N);
% figure('Position',[2000 10 600 600])
% hold on
% plot(xplot,F_sp,'linewidth',2)
% plot(xplot, Fc,'--','color',[.6 .6 .6],'linewidth',2)
% plot(xplot,-Fs,'--','color',[.8 .8 .8],'linewidth',2)
% legend('spatial filter','center','surround','Location','southeast')
% xlabel('X (deg)')
% title([info.cell.type, ' ', info.cell.rftype, ' cell'])
% set(gca,'fontsize',18,'linewidth',2)
% box off