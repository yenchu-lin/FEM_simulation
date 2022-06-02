function [info] = femc_eye_traj3(info)
% Load eye trajectory
% Modified for ddpi data 
% Last updated on 2021/09/21 by YCL

% load the single eye trajectory
datadir = info.eyemvnt.datadir;
datafile = [datadir filesep info.eyemvnt.filename];
eyetraj = importdata(datafile,'eyetraj');
fname = info.eyemvnt.fname;
info.eyemvnt.ntraj = length(eyetraj.(fname));
conversionFactor = 1000 / eyetraj.(fname){1}.sRate;

% define the number of timepoints to use from the eye trajectory
duration = info.eyemvnt.t.ramp + info.eyemvnt.t.plateau + info.eyemvnt.t.post;
info.eyemvnt.t.duration = duration;
info.eyemvnt.nsamp = floor(duration/conversionFactor);

% center of the image
new_center = (info.img.imsize*info.img.npad)/2;

% get a list of trials to be included (some trials end early)
ilist = [];
for itraj = 1:info.eyemvnt.ntraj
    
    % save the eye positions from start to end
    tstart = round(eyetraj.(fname){itraj}.TimeHoldOFF / conversionFactor);
    tend = tstart + info.eyemvnt.nsamp - 1;
    X = [];
    if size(eyetraj.(fname){itraj}.xshift,2) >=  tend
        ilist = [ilist,itraj];
    end
end

emv.x       = zeros(info.eyemvnt.nsamp,2,length(ilist),info.cell.ncell);
emv.xrotate = zeros(info.eyemvnt.nsamp,2,length(ilist),info.cell.ncell);
for ii = 1:length(ilist)
    itraj = ilist(ii);
    
    tstart = round(eyetraj.(fname){itraj}.TimeHoldOFF / conversionFactor);
    tend = tstart + info.eyemvnt.nsamp - 1;
    x = eyetraj.(fname){itraj}.xshift(tstart:tend); % in arcmin
    y = eyetraj.(fname){itraj}.yshift(tstart:tend); % in arcmin

%     % s-golay filtering 
%     smoothing = 31;
%     smx = sgfilt(x, 3, smoothing, 0);
%     smy = sgfilt(y, 3, smoothing, 0);
% 
%     X(:,1) = smx;
%     X(:,2) = smy;
    
    X(:,1) = x;
    X(:,2) = y;

    % rescale and recenter eye trajectory 
    % eye movement from the data is viewing 2 deg images
    pixel_deg = info.img.imsize / info.img.D; % number of pixel/deg
    pixel_arcmin = pixel_deg / 60; % number of pixel/arcmin
    scaledX = X * pixel_arcmin + new_center; % X is in arcmin; scaledX is in pixel 

    for icell = 1:info.cell.ncell
        newX = scaledX + repmat(info.cell.ctr(icell,:)*pixel_arcmin,length(X),1);    
        if info.eyemvnt.rotate == 0
            emv.x(:,:,ii,icell) = newX;
        else
            xrotate = imrotate(newX, info.eyemvnt.rotate);
            emv.x(:,:,ii,icell) = xrotate';
        end
    end

end

info.emv = emv;
info.eyemvnt.ntraj = length(ilist);
info.eyemvnt.conversionFactor = conversionFactor;

end

