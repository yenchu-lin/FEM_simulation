function [info] = femc_cell_center(info,model)
% Hexagonal array of cells based on the radius
% The function sets the centers of the cells based on the radius of the
% center RF of each cell type

% How many times of radius we set it up as the distance between cells
scale_p = 4;
scale_m = 2;

% Set the centers of the cells based on r (radius of the center RF of the cell)
if strcmpi(info.cell.type,'p')
    info.cell.spat.celltype_no = 1; % 1- P cell
    % The cell array has 12 columns and 15 rows = 180 cells
    ncol = 12;
    nrow = 15;
    r = model{1}.component{info.cell.spat.celltype_no}.rc / sqrt(2) * 60 * scale_p; % in arcmin
elseif strcmpi(info.cell.type,'m')
    info.cell.spat.celltype_no = 2; % 2- M cell
    % The cell array has 8 columns and 10 rows = 80 cells
    ncol = 8;
    nrow = 10;
    r = model{1}.component{info.cell.spat.celltype_no}.rc / sqrt(2) * 60 * scale_m; % in arcmin
end

n = ncol*nrow;
w = sqrt(3)*r; % width of the hexagon; also as horizontal spacing
h = 2*r;       % height of the hexagon
v = h * 3/4;   % vertical spacing

% irregularity of the array
mu = 0; % arcmin
sd = 1.3;
dx = mu + sd .* randn(n,1);
dy = mu + sd .* randn(n,1);

ctr = zeros(n,2);
ctr2 = zeros(n,2);
for ir = 1:nrow
    for ic = 1:ncol
        icell = (ir-1)*ncol + ic;
        if mod(ir,2)
            ctr(icell,:)  = [(ic-1)*w,(ir-1)*v];
            ctr2(icell,:) = [(ic-1)*w+dx(icell),(ir-1)*v+dy(icell)];
        else
            ctr(icell,:)  = [(ic-1)*w+(w/2),(ir-1)*v];
            ctr2(icell,:) = [(ic-1)*w+(w/2)+dx(icell),(ir-1)*v+dy(icell)];
        end
    end
end

total_width  = w*(ncol-1)+w/2;
total_height = v*(nrow-1);
new_center = [-total_width/2, -total_height/2];

% use the array with regularity
ctr = ctr + repmat(new_center, n, 1);
% randcell = randi(n, info.cell.ncell, 1); % with repeats
randcell = randperm(n, info.cell.ncell); % without repeats
info.cell.ctr = ctr(randcell,:);

% use the array with irregularity
% ctr2 = ctr2 + repmat(new_center,n,1);
% info.cell.ctr = ctr2(randcell,:);
% info.cell.ctr_array = ctr(randcell,:);

end