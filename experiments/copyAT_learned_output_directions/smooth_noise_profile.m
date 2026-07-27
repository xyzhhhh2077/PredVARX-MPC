function [sigma_w, sigma_e] = smooth_noise_profile(T, cycle, sw_min, sw_max, se_min, se_max, phase_e)
%SMOOTH_NOISE_PROFILE Raised-cosine standard-deviation envelopes.
% Noise samples remain Gaussian white; only their instantaneous scale varies.
validateattributes(T, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cycle, {'numeric'}, {'scalar','integer','positive'});
assert(sw_min >= 0 && sw_max >= sw_min, 'Invalid process-noise range.');
assert(se_min >= 0 && se_max >= se_min, 'Invalid measurement-noise range.');

k = 0:T-1;
theta = 2*pi*k/cycle;
sigma_w = sw_min + 0.5*(sw_max-sw_min).*(1-cos(theta));
sigma_e = se_min + 0.5*(se_max-se_min).*(1-cos(theta+phase_e));
end
