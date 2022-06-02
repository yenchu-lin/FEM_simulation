function [p] = femc_lgn_multiplier(RF,model,info,filtertype)

model_no = info.cell.(filtertype).model_no;
celltype_no = info.cell.(filtertype).celltype_no;
rftype = info.cell.rftype;
celltype = info.cell.type;

% Spatial type 1: circular Gaussian
if model_no == 1
    if strcmpi(info.cell.type,'p')
        label_no = 1;
    elseif strcmpi(info.cell.type,'m')   
        label_no = 2;
    end
    if strcmp(RF, 'center') 
        Kc       = model{model_no}.component{label_no}.Kc;
        rc       = model{model_no}.component{label_no}.rc;
        p.str    = Kc * pi * (rc^2); % A; strength
        p.radius = rc/(sqrt(2)); % sigma (in deg); radius
        p.Kc = Kc;
        p.rc = rc;
    elseif strcmp(RF, 'surround')
        Ks       = model{model_no}.component{label_no}.Ks;
        rs       = model{model_no}.component{label_no}.rs;
        p.str    = Ks * pi * (rs^2);
        p.radius = rs/(sqrt(2)); % in deg
        p.Ks = Ks;
        p.rs = rs;
    end
 
% Temporal type 1: gamma difference
elseif model_no == 2 
    p.k      = model{model_no}.component{celltype_no}.k;
    p.c      = model{model_no}.component{celltype_no}.c;
    p.t0     = model{model_no}.component{celltype_no}.t0;
    p.n      = model{model_no}.component{celltype_no}.n;
    p.tdelay = model{model_no}.component{celltype_no}.tdelay;

% Temporal type 2: HPLP
elseif model_no == 3
    if strcmpi(celltype,'P') && strcmpi(rftype,'ON')
        if strcmp(RF, 'center')
            label_no = 1;
        elseif strcmp(RF, 'surround')
            label_no = 2;
        end
    elseif strcmpi(celltype,'P') && strcmpi(rftype,'OFF')
        if strcmp(RF, 'center')
            label_no = 3;
        elseif strcmp(RF, 'surround')
            label_no = 4;
        end
    elseif strcmpi(celltype,'M') && strcmpi(rftype,'ON')
        label_no = 5;
    elseif  strcmpi(celltype,'M') && strcmpi(rftype,'OFF')
        label_no = 6;
    end
    p.D      = model{model_no}.component{label_no}.D;
    p.A      = model{model_no}.component{label_no}.A;
    p.Hs     = model{model_no}.component{label_no}.Hs;
    p.tau_S  = model{model_no}.component{label_no}.tau_S;
    p.tau_L  = model{model_no}.component{label_no}.tau_L;
    p.NL     = model{model_no}.component{label_no}.NL;
    p.tdelay = model{model_no}.component{label_no}.tdelay;
end

end
