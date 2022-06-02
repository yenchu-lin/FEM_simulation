function [model,info] = femc_lgn_parameters(info)
% The function returns the parameters for different LGN models

if strcmpi(info.cell.type,'p')
    info.cell.spat.celltype_no = 1; % 1 - P cell; 2 - M cell
    if info.cell.temp.model_no == 2
        info.cell.temp.celltype_no = 1; % 1- P cell; 2- M cell
    elseif info.cell.temp.model_no == 3
        info.cell.temp.celltype_no = 1; % 1- P ON cell or 3- P OFF cell
    end
elseif strcmpi(info.cell.type,'m')
    info.cell.spat.celltype_no = 2;
    if info.cell.temp.model_no == 2
        info.cell.temp.celltype_no = 2; % 1- P cell; 2- M cell
    elseif info.cell.temp.model_no == 3
        info.cell.temp.celltype_no = 5; % 5- M ON cell or 6- M cell OFF cell
    end
end

% P and M cells; On and off cells in primate
% -------------------------------------------------
% Spatial 
% type1: circular Gaussian 
%       (Ref: Croner and Kaplan (1994) Receptive fields of P and M ganglion
%        cells across the primate retina)
% -------------------------------------------------
% Temporal
% type1: Gamma function model (Ref: Rucci 2000)
% type2: HPLP model 
%        (Ref: Victor 1987 & Benardete and Kaplan 1999; M On cell in primates & Benardete and Kaplan 1997 P cell in primates)
% -------------------------------------------------

% Spatial type 1: circular Gaussian
% Ref: Croner and Kaplan (1995) Receptive fields of P and M ganglion
%      cells across the primate retina
model{1}.name = 'spatial_circular_Gaussian';
model{2}.name = 'temporal_gamma';
model{3}.name = 'temporal_hplp';

model{1}.component{1}.type = 'P cell';
model{1}.component{1}.rc = 0.03;
model{1}.component{1}.Kc = 325.2;
model{1}.component{1}.rs = 0.18;
model{1}.component{1}.Ks = 4.4;

model{1}.component{2}.type = 'M cell';
model{1}.component{2}.rc = 0.10;
model{1}.component{2}.Kc = 148.0;
model{1}.component{2}.rs = 0.72;
model{1}.component{2}.Ks = 1.1;


% Temporal type 1: gamma (need to be corrected)
model{2}.component{1}.type = 'P cell';
model{2}.component{1}.k = [1;0.6];
model{2}.component{1}.c = [60;40]; % 1/sec
model{2}.component{1}.t0 = [0;0];
model{2}.component{1}.n = [2;2];
model{2}.component{1}.tdelay = 3/1000; %sec

model{2}.component{2}.type = 'M cell';
model{2}.component{2}.k = [1;0.6];
model{2}.component{2}.c = [60;40]; % 1/sec
model{2}.component{2}.t0 = [0;0];
model{2}.component{2}.n = [2;2];
model{2}.component{2}.tdelay = 3/1000; %sec

% Temporal type 2: HPLP
cnumber = 1;
model{3}.component{cnumber}.type = 'P ON cell center';
model{3}.component{cnumber}.D = 3.5/1000; % an initial delay
model{3}.component{cnumber}.A = 67.59;% the overall gain
model{3}.component{cnumber}.Hs = 0.69; % the strength of the high-pass stage
model{3}.component{cnumber}.tau_S = 29.36/1000; % the time constant of the high-pass stages; tau0
model{3}.component{cnumber}.tau_L = 1.27/1000; % the time constant of the low-pass stages
model{3}.component{cnumber}.NL = 38; % the number of low-pass stages
model{3}.component{cnumber}.tdelay = 3/1000; % sec

cnumber = 2;
model{3}.component{cnumber}.type = 'P ON cell surround';
model{3}.component{cnumber}.D = 3.5/1000; % an initial delay
model{3}.component{cnumber}.A = 49.98;% the overall gain
model{3}.component{cnumber}.Hs = 0.48; % the strength of the high-pass stage
model{3}.component{cnumber}.tau_S = 18.62/1000; % the time constant of the high-pass stages; tau0
model{3}.component{cnumber}.tau_L = 0.50/1000; % the time constant of the low-pass stages
model{3}.component{cnumber}.NL = 111; % the number of low-pass stages
model{3}.component{cnumber}.tdelay = 7.24/1000; % sec

cnumber = 3;
model{3}.component{cnumber}.type = 'P OFF cell center';
model{3}.component{cnumber}.D = 3.5/1000; % an initial delay
model{3}.component{cnumber}.A = 52.95;% the overall gain
model{3}.component{cnumber}.Hs = 0.75; % the strength of the high-pass stage
model{3}.component{cnumber}.tau_S = 35.68/1000; % the time constant of the high-pass stages; tau0
model{3}.component{cnumber}.tau_L = 1.87/1000; % the time constant of the low-pass stages
model{3}.component{cnumber}.NL = 27; % the number of low-pass stages
model{3}.component{cnumber}.tdelay = 3/1000; % sec

cnumber = 4;
model{3}.component{cnumber}.type = 'P OFF cell surround';
model{3}.component{cnumber}.D = 3.5/1000; % an initial delay
model{3}.component{cnumber}.A = 23.91;% the overall gain
model{3}.component{cnumber}.Hs = 0.42; % the strength of the high-pass stage
model{3}.component{cnumber}.tau_S = 52.00/1000; % the time constant of the high-pass stages; tau0
model{3}.component{cnumber}.tau_L = 0.90/1000; % the time constant of the low-pass stages
model{3}.component{cnumber}.NL = 67; % the number of low-pass stages
model{3}.component{cnumber}.tdelay = 10.22/1000; % sec

cnumber = 5;
model{3}.component{cnumber}.type = 'M ON cell';
model{3}.component{cnumber}.D = 2/1000; % an initial delay
model{3}.component{cnumber}.A = 499.77;% the overall gain
model{3}.component{cnumber}.Hs = 1.00; % the strength of the high-pass stage
model{3}.component{cnumber}.tau_S = 37.34/1000; % the time constant of the high-pass stages; tau0
model{3}.component{cnumber}.tau_L = 1.34/1000; % the time constant of the low-pass stages
model{3}.component{cnumber}.NL = 30; % the number of low-pass stages
model{3}.component{cnumber}.tdelay = 3/1000; % sec

cnumber = 6;
model{3}.component{cnumber}.type = 'M OFF cell';
model{3}.component{cnumber}.D = 2.4/1000; % an initial delay
model{3}.component{cnumber}.A = 517.40;% the overall gain
model{3}.component{cnumber}.Hs = 0.97; % the strength of the high-pass stage
model{3}.component{cnumber}.tau_S = 57.62/1000; % the time constant of the high-pass stages; tau0
model{3}.component{cnumber}.tau_L = 1.97/1000; % the time constant of the low-pass stages
model{3}.component{cnumber}.NL = 22; % the number of low-pass stages
model{3}.component{cnumber}.tdelay = 3/1000; % sec
