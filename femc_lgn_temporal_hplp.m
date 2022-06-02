function [X,K,t] = femc_lgn_temporal_hplp(p, delta_t, n_half)

% temporal component of LGN model
  % high-pass low-pass model
  % Ref: Victor 1987 & Benardete and Kaplan 1997; 1999
  % input: p -- a struture of the parameters
  %        delta_T -- the time step (inverse of the sampling rate)
  %        n_half -- the step size for omega
  % return: the time series

fmax = 1/2/delta_t; % Nyquist frequency (half the sampling rate); in Hz
wmax = fmax*2*pi; % Nyquist frequency in radian
wmin = wmax/n_half; % in radian
w = (0:n_half)*wmin; % the list of w; in radian

% in freq domain
K_list = p.A * exp(-1i * w * p.D) .* (1 - p.Hs ./ (1 + 1i * w * p.tau_S))...
    .* (1 ./ (1 + 1i * w * p.tau_L).^ p.NL); % radian/sec

% % plot
% figure
% loglog(w/2/pi, abs(K_list)) % plot the magnitude of K
% ylim([p.A*0.001 p.A])
% xlabel('Frequency (Hz)')
% ylabel('Response (impulses/sec)')

K = [K_list(1:end-1),real(K_list(end)),fliplr(conj(K_list(2:end-1)))]; % K(w)
X = ifft(K) * wmin * 2 * n_half / (2*pi); 
t = (0:length(K)-1) .* delta_t;

% transfer function
% N = p.NL - 1;
% H = t.^(N) ./ (factorial(N) * (p.tau_L^p.NL)) .* exp(-t/p.NL);

end