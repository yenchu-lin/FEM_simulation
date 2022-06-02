function [info] = femc_cell_center_old(info,model)
% [imginfo] = femc_cell_center_old(imginfo,model,smodel_no,scelltype_no)

% Set the centers of the cells based on r (radius of the center RF of the cell)
if strcmpi(info.cell.type,'p')
    info.cell.spat.celltype_no = 1; % 1- P cell
    r = model{1}.component{info.cell.spat.celltype_no}.rc / sqrt(2) * 60 * 4; % in arcmin
elseif strcmpi(info.cell.type,'m')
    info.cell.spat.celltype_no = 2; % 2- M cell
    r = model{1}.component{info.cell.spat.celltype_no}.rc / sqrt(2) * 60 * 4; % in arcmin
end

% % set the center of the cells
% r = model{smodel_no}.component{scelltype_no}.rc / sqrt(2) * 60; % in arcmin

% hexagonal based on the radius
imginfo.ctr(1,:) = [0,0];
imginfo.ctr(2,:) = [0,sqrt(3)*r];
imginfo.ctr(3,:) = [0,-sqrt(3)*r];
imginfo.ctr(4,:) = [1.5*r,sqrt(3)/2*r];
imginfo.ctr(5,:) = [1.5*r,-sqrt(3)/2*r];
imginfo.ctr(6,:) = [-1.5*r,sqrt(3)/2*r];
imginfo.ctr(7,:) = [-1.5*r,-sqrt(3)/2*r];
imginfo.ctr(8,:) = [3*r,0];
imginfo.ctr(9,:) = [-3*r,0];
imginfo.ctr(10,:) = [3*r,sqrt(3)*r];
imginfo.ctr(11,:) = [3*r,-sqrt(3)*r];
imginfo.ctr(12,:) = [-3*r,sqrt(3)*r];
imginfo.ctr(13,:) = [-3*r,-sqrt(3)*r];
imginfo.ctr(14,:) = [1.5*r,3*sqrt(3)/2*r];
imginfo.ctr(15,:) = [1.5*r,-3*sqrt(3)/2*r];
imginfo.ctr(16,:) = [-1.5*r,3*sqrt(3)/2*r];
imginfo.ctr(17,:) = [-1.5*r,-3*sqrt(3)/2*r];
imginfo.ctr(18,:) = [0,2*sqrt(3)*r];
imginfo.ctr(19,:) = [0,-2*sqrt(3)*r];
imginfo.ctr(20,:) = [3*r,2*sqrt(3)*r];
imginfo.ctr(21,:) = [3*r,-2*sqrt(3)*r];
imginfo.ctr(22,:) = [-3*r,2*sqrt(3)*r];
imginfo.ctr(23,:) = [-3*r,-2*sqrt(3)*r];
imginfo.ctr(24,:) = [1.5*r,5*sqrt(3)/2*r];
imginfo.ctr(25,:) = [1.5*r,-5*sqrt(3)/2*r];
imginfo.ctr(26,:) = [-1.5*r,5*sqrt(3)/2*r];
imginfo.ctr(27,:) = [-1.5*r,-5*sqrt(3)/2*r];
imginfo.ctr(28,:) = [0,3*sqrt(3)*r];
imginfo.ctr(29,:) = [0,-3*sqrt(3)*r];
imginfo.ctr(30,:) = [4.5*r,sqrt(3)/2*r];
imginfo.ctr(31,:) = [4.5*r,-sqrt(3)/2*r];
imginfo.ctr(32,:) = [-4.5*r,sqrt(3)/2*r];
imginfo.ctr(33,:) = [-4.5*r,-sqrt(3)/2*r];
imginfo.ctr(34,:) = [4.5*r,sqrt(3)*3/2*r];
imginfo.ctr(35,:) = [4.5*r,-sqrt(3)*3/2*r];
imginfo.ctr(36,:) = [-4.5*r,sqrt(3)*3/2*r];
imginfo.ctr(37,:) = [-4.5*r,-sqrt(3)*3/2*r];
imginfo.ctr(38,:) = [4.5*r,sqrt(3)*5/2*r];
imginfo.ctr(39,:) = [4.5*r,-sqrt(3)*5/2*r];
imginfo.ctr(40,:) = [-4.5*r,sqrt(3)*5/2*r];
imginfo.ctr(41,:) = [-4.5*r,-sqrt(3)*5/2*r];
imginfo.ctr(42,:) = [5.5*r,0];
imginfo.ctr(43,:) = [-5.5*r,0];
imginfo.ctr(44,:) = [5.5*r,sqrt(3)*r];
imginfo.ctr(45,:) = [5.5*r,-sqrt(3)*r];
imginfo.ctr(46,:) = [-5.5*r,sqrt(3)*r];
imginfo.ctr(47,:) = [-5.5*r,-sqrt(3)*r];
imginfo.ctr(48,:) = [5.5*r,2*sqrt(3)*r];
imginfo.ctr(49,:) = [5.5*r,-2*sqrt(3)*r];
imginfo.ctr(50,:) = [-5.5*r,2*sqrt(3)*r];
imginfo.ctr(51,:) = [-5.5*r,-2*sqrt(3)*r];
imginfo.ctr(52,:) = [3*r,3*sqrt(3)*r];
imginfo.ctr(53,:) = [3*r,-3*sqrt(3)*r];
imginfo.ctr(54,:) = [-3*r,3*sqrt(3)*r];
imginfo.ctr(55,:) = [-3*r,-3*sqrt(3)*r];
imginfo.ctr(56,:) = [1.5*r,7*sqrt(3)/2*r];
imginfo.ctr(57,:) = [1.5*r,-7*sqrt(3)/2*r];
imginfo.ctr(58,:) = [-1.5*r,7*sqrt(3)/2*r];
imginfo.ctr(59,:) = [-1.5*r,-7*sqrt(3)/2*r];
imginfo.ctr(60,:) = [4.5*r,sqrt(3)*7/2*r];
imginfo.ctr(61,:) = [4.5*r,-sqrt(3)*7/2*r];
imginfo.ctr(62,:) = [-4.5*r,sqrt(3)*7/2*r];
imginfo.ctr(63,:) = [-4.5*r,-sqrt(3)*7/2*r];

info.cell.ctr = imginfo.ctr;

% random
% for icell = 1:imginfo.ncell
%     if icell == 1
%         ctr = [0,0];
%     else
%         ctr = [randi([-30 30],1),randi([-30 30],1)]; % in arcmin
%     end
%    imginfo.ctr(icell,:) = ctr;
% end

end

