function [ts_ramp,ramp_contrast] = femc_ramp(ts_spatial,info)
% ramp protocol for the contrast

Tup = 1;
Tsteady = Tup + round(info.eyemvnt.t.ramp / info.eyemvnt.conversionFactor);
Tdown = Tsteady + round(info.eyemvnt.t.plateau / info.eyemvnt.conversionFactor);
Toff = round(info.eyemvnt.t.duration / info.eyemvnt.conversionFactor);

ramp_contrast = zeros(1,Toff);
ramp_contrast(Tup:Tsteady-1) = linspace(0,1,Tsteady-Tup);
ramp_contrast(Tsteady:Tdown-1) = 1;
ramp_contrast(Tdown:Toff-1) = linspace(1,0,Toff-Tdown);
ramp_contrast = ramp_contrast';

ts_ramp = ts_spatial .* repmat(ramp_contrast,1,2,info.img.nstim);

end

