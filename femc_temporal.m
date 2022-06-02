function [Timeseries,info] = femc_temporal(model,ts_ramp,info)

% Cell type
% HPLP model: 1- P ON cell or 3- P OFF cell or 5- M ON cell or 6 - M cell OFF cell
% Gamma model: 1- P cell, 2 - M cell

% if strcmpi(info.cell.type,'p')
%     if info.cell.temp.model_no == 2
%         info.cell.temp.celltype_no = 1;
%     elseif info.cell.temp.model_no == 3
%         if strcmp(info.cell.rftype,'ON') == 1
%             info.cell.temp.celltype_no = 1;
%         elseif strcmp(info.cell.rftype,'OFF') == 1
%             info.cell.temp.celltype_no = 3;
%         end
%     end
%  
% elseif strcmpi(info.cell.type,'m')
%     if info.cell.temp.model_no == 2
%         info.cell.temp.celltype_no = 2; 
%     elseif info.cell.temp.model_no == 3
%         if strcmp(info.cell.rftype,'ON') == 1
%             info.cell.temp.celltype_no = 5;
%         elseif strcmp(info.cell.rftype,'OFF') == 1
%             info.cell.temp.celltype_no = 6;
%         end
%     end
% end

% Temporal filtering
nt = length(ts_ramp);
P = nextpow2(nt);
nsamp = 2^(P+1); % power of 2 (for the benefit of computational power)
Timeseries = zeros(nt,info.img.nstim);
for istim = 1:info.img.nstim
    RF = 'center';
    [Gcent,~] = femc_temporal_filter(ts_ramp(:,1,istim), RF, model, nsamp, info);
    RF = 'surround';
    [Gsurr,~] = femc_temporal_filter(ts_ramp(:,2,istim), RF, model, nsamp, info);
    G = Gcent - Gsurr;
    G(G<0) = 0;
    Timeseries(:,istim) = G(1:nt);
end

end
