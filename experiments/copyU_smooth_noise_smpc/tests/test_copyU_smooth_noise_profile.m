%% test_copyU_smooth_noise_profile
% RED/GREEN behavior test for a smooth nonstationary Gaussian noise envelope.
clear; clc; addpath(fileparts(fileparts(mfilename('fullpath'))));

T = 1200;
cycle = 400;
[sigma_w, sigma_e] = smooth_noise_profile(T, cycle, 0.02, 0.09, 0.025, 0.10, pi/3);

assert(isrow(sigma_w) && numel(sigma_w) == T, 'sigma_w must be 1-by-T');
assert(isrow(sigma_e) && numel(sigma_e) == T, 'sigma_e must be 1-by-T');
assert(abs(min(sigma_w) - 0.02) < 1e-12, 'sigma_w minimum is wrong');
assert(abs(max(sigma_w) - 0.09) < 1e-12, 'sigma_w maximum is wrong');
assert(abs(min(sigma_e) - 0.025) < 2e-6, 'sigma_e minimum is wrong');
assert(abs(max(sigma_e) - 0.10) < 2e-6, 'sigma_e maximum is wrong');
assert(max(abs(diff(sigma_w,2))) < 2e-5, 'sigma_w envelope is not smooth enough');
assert(max(abs(diff(sigma_e,2))) < 2e-5, 'sigma_e envelope is not smooth enough');
assert(norm(sigma_w(1:cycle)-sigma_w(cycle+1:2*cycle),inf) < 1e-12, 'sigma_w is not periodic');
assert(norm(sigma_e(1:cycle)-sigma_e(cycle+1:2*cycle),inf) < 1e-12, 'sigma_e is not periodic');
assert(std(sigma_w) > 0.02 && std(sigma_e) > 0.02, 'envelopes must vary materially');

fprintf('PASS test_copyU_smooth_noise_profile: sw=[%.3f,%.3f], se=[%.3f,%.3f]\n', ...
    min(sigma_w), max(sigma_w), min(sigma_e), max(sigma_e));
