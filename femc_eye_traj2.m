function [info] = femc_eye_traj2(info)
% Load eye trajectory

% load the single eye trajectory
datadir = info.eyemvnt.datadir;
datafile = [datadir filesep info.eyemvnt.filename];
eyetraj = importdata(datafile,'eyetraj');
fname = info.eyemvnt.fname;
info.eyemvnt.ntraj = length(eyetraj.(fname));

% define the number of timepoints to use from the eye trajectory
duration = info.eyemvnt.t.ramp + info.eyemvnt.t.plateau + info.eyemvnt.t.post;
info.eyemvnt.t.duration = duration;
info.eyemvnt.nsamp = duration;

% center of the image
new_center = (info.img.imsize*info.img.npad)/2;

emv.x = zeros(info.eyemvnt.nsamp,2,info.eyemvnt.ntraj,info.cell.ncell);
emv.xrotate = zeros(info.eyemvnt.nsamp,2,info.eyemvnt.ntraj,info.cell.ncell);
for itraj = 1:info.eyemvnt.ntraj
    
    % save the eye positions from start to end
    tstart = round(eyetraj.(fname){itraj}.TimeHoldOFF);
    tend = tstart + info.eyemvnt.nsamp - 1;
    X = [];
    x = eyetraj.(fname){itraj}.xshift(tstart:tend); % in arcmin
    y = eyetraj.(fname){itraj}.yshift(tstart:tend); % in arcmin
    
    % s-golay filtering 
    smoothing = 31;
    smx = sgfilt(x, 3, smoothing, 0);
    smy = sgfilt(y, 3, smoothing, 0);
    
    X(:,1) = smx;
    X(:,2) = smy;
    
    % rescale and recenter eye trajectory 
    % eye movement from the data is viewing 2 deg images
    pixel_deg = info.img.imsize / info.img.D; % number of pixel/deg
    pixel_arcmin = pixel_deg / 60;
    scaledX = X * pixel_arcmin + new_center;
    
    for icell = 1:info.cell.ncell
        newX = scaledX + repmat(info.cell.ctr(icell,:)*pixel_arcmin,length(X),1);    
        if info.eyemvnt.rotate == 0
            emv.x(:,:,itraj,icell) = newX;
        else
            xrotate = imrotate(newX, info.eyemvnt.rotate);
            emv.x(:,:,itraj,icell) = xrotate';
        end
    end
end

info.emv = emv;

% figure
% subplot(1,3,1)
% plot(X(:,1), X(:,2))
% xlabel('x (deg)')
% ylabel('y (deg)')
% axis square
% axis equal
% subplot(1,3,2)
% plot(scaledX(:,1),scaledX(:,2))
% xlabel('x (pixel)')
% ylabel('y (pixel)')
% axis square
% axis equal
% title('rescaled')
% subplot(1,3,3)
% plot(newX(:,1),newX(:,2))
% xlabel('x (pixel)')
% ylabel('y (pixel)')
% axis square
% axis equal
% title('recentered')

end

