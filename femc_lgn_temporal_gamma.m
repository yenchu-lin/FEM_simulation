function [G] = femc_lgn_temporal_gamma(p,t)
% temporal component of LGN model
% difference of two gamma functions
% Ref: Rucci 2000

g = cell(2,1);
for ii=1:2
    g{ii} = (((p.c(ii)*(t-p.t0(ii))).^p.n(ii)) .* exp(-p.c(ii)*(t-p.t0(ii))))/...
        (p.n(ii)^p.n(ii)*exp(-p.n(ii)));
end

G = p.k(1)*g{1} - p.k(2)*g{2};

end