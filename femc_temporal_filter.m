function [G,K,t] = femc_temporal_filter(time_series, RF, model, nsamp, info)

% Add npad timepoints
n = length(time_series);
npad = nsamp - n; 
input_time_series = [time_series',zeros(1,npad)];  % power of 2 (for the benefit of computational power)
n_half = length(input_time_series)/2;
FT = fft(input_time_series);

filtertype = 'temp';
[p] = femc_lgn_multiplier(RF,model,info,filtertype);

% HPLP Model
if info.cell.temp.model_no == 3 
    if strcmp(RF,'center')
        [~,K,~] = femc_lgn_temporal_hplp(p, info.cell.temp.delta_t, n_half);
        F = K .* FT;
        G = ifft(F); % back to the time domain
        t = (0:length(G)-1) .* info.cell.temp.delta_t;
    elseif strcmp(RF,'surround')
        if info.cell.temp.celltype_no == 1 || info.cell.temp.celltype_no == 3  % P cell
            info.cell.temp.celltype_no = info.cell.temp.celltype_no + 1;
            [p] = femc_lgn_multiplier(RF,model,info,filtertype);
            [~,K,~] = femc_lgn_temporal_hplp(p, info.cell.temp.delta_t, n_half);
            F = K .* FT;
            G = ifft(F); % back to the time domain
            t = (0:length(G)-1) .* info.cell.temp.delta_t;
        else % M cell
            [~,K,~] = femc_lgn_temporal_hplp(p, info.cell.temp.delta_t, n_half);
            F = K .* FT;
            G = ifft(F);
            t = (0:length(G)-1) .* info.cell.temp.delta_t;
            G = G .* (t-p.tdelay);
        end
    end

% Difference of Two Gamma Model
elseif info.cell.temp.model_no == 2
    t_list = 0:delta_t:n_half*delta_t;
    [K] = femc_lgn_temporal_gamma(p,t_list);
    K_list = fft(K); 
    FT_K = [K_list(1:end-1),real(K_list(end)),fliplr(conj(K_list(2:end-1)))];
    F = FT_K .* FT;
    G = ifft(F);  
    t = (0:length(G)-1) .* info.cell.temp.delta_t;   
    if strcmp(RF,'surround')
        G = G .* (t-p.tdelay);
    end
end

end

