function [info] = femc_input_eye_movement(info)

% obtain the eye movements from experimental data
if strcmp(info.eyemvnt.type,'exp_traj')
    if ~info.eyemvnt.ddpi
        [info] = femc_eye_traj2(info);
    elseif info.eyemvnt.ddpi
        [info] = femc_eye_traj3(info);
    end

% adding a simulated trajectory for testing purpose (only two horizontal traj   
elseif strcmp(info.eyemvnt.type,'sim_traj')
    
    % image size = 512 x 512 pixels without zero padding (npad = 2) 
    % grating range from 256 - 768
    lengthx = 1500;
    info.eyemvnt.ntraj = 3;
    info.eyemvnt.conversionFactor = 1;
    info.eyemvnt.t.duration = info.eyemvnt.t.ramp + info.eyemvnt.t.plateau + info.eyemvnt.t.post;
    for icell = 1:info.cell.ncell
        for itraj = 1:info.eyemvnt.ntraj
            info.emv.x(:,:,itraj,icell) = ...
                [linspace(300,720,1500);repmat(500+50*(icell-1),1,lengthx)]' ...
                + repmat(info.cell.ctr(icell,:),lengthx,1);
        end
    end

% obtain the eye movements from simulated data
% elseif strcmp(info.eyemvnt.type,'sim_traj') == 1
%     
%     % more input
%     info.eyemvnt.filename = '';
%     info.eyemvnt.ntraj = 10;
%     info.eyemvnt.T = 1;
%     info.eyemvnt.cov = [0.8,0.1;0.1,0.2];
%     info.eyemvnt.dt = 0.001;
%     info.eyemvnt.nexamp = 1000;
%     
%     [emv,info.img] = femc_eye_traj_sim(info.eyemvnt, info.img);
%     
% end

end
