function out = evaluate_multistep_prediction(y_raw,u_raw,run_id,model,scales,H,start_step)
%EVALUATE_MULTISTEP_PREDICTION Recursive multi-step prediction on real data.
% Benchmark protocol mirrors Mohajerin & Waslander ICRA 2018 (arXiv
% 1806.00526): fixed init window, then iterate the learned model H steps
% with the recorded input sequence; report per-horizon error distributions
% over all start points (paper reports mean + 99th percentile).
%
% Units (Pelican, verified): Pos m, Euler rad, Vel m/s = 100*diff(Pos),
% pqr rad/s (body frame). Vel/pqr are DERIVED so they are not in y; we
% reconstruct speed from predicted positions the same way the dataset does.
%
% Inputs:
%   y_raw     p x T  raw measurements (Pos m, Euler rad, Motors)
%   u_raw     m x T  raw commands
%   run_id    1 x T  segment ids (transitions only inside a segment)
%   model     struct with A,B,P,R,y_mean,u_mean (trained in STANDARDIZED space)
%   scales    struct y_offset,y_scale,u_offset,u_scale (training statistics)
%   H         prediction horizon in steps (paper: 40 = 0.4s, 190 = 1.9s)
%   start_step sampling stride for start points (default H, non-overlap)
%
% Outputs:
%   err_pos  [N x H]   position error norm  (m)   -> cm in caller
%   err_att  [N x H]   Euler angle error norm (rad)
%   err_vel  [N x H-1] speed error norm (m/s):
%                     |100*(Pos_pred(h+1)-Pos_pred(h)) - 100*(Pos_raw(k+h+1)-Pos_raw(k+h))|
%   err_rate [N x H-1] Euler-rate error norm (rad/s)
%   err_rate_body [N x H-1] body-frame rate error via true-attitude transform
%                     (closest fair comparison to the paper's pqr metric)
%   per_*    same quantities for persistence (hold last observation)
%   starts  1 x N  start indices used
%   N        number of starts
if nargin < 7 || isempty(start_step), start_step = H; end
[p,T] = size(y_raw); m = size(u_raw,1);
run_id = run_id(:)';
assert(size(u_raw,2)==T && numel(run_id)==T,'Sample counts must match.');
assert(p==10,'Expected Pelican 10-channel raw y.');
H = min(H,T-2);

% Standardize with TRAINING statistics only (same as the fit).
yc = (y_raw - scales.y_offset)./scales.y_scale;
uc = (u_raw - scales.u_offset)./scales.u_scale;

% Candidate start points: inside a segment, at least 1 sample of history,
% prediction window [k0+1, k0+H] must stay inside the segment; spacing
% start_step (default H, i.e. non-overlapping windows).
changes = find(diff(run_id)~=0);            % last index of each segment
seg_start = [1, changes+1];
seg_end   = [changes, T];
starts = [];
for s = 1:numel(seg_start)
    a = seg_start(s); b = seg_end(s);
    for k0 = a+1 : start_step : b-H
        starts(end+1) = k0; %#ok<AGROW>
    end
end
N = numel(starts);
assert(N>0,'No valid start points for this segment layout.');

err_pos  = zeros(N,H); err_att  = zeros(N,H);
err_vel  = zeros(N,H-1); err_rate = zeros(N,H-1); err_rate_body = zeros(N,H-1);
per_pos  = zeros(N,H); per_att  = zeros(N,H);
per_vel  = zeros(N,H-1); per_rate = zeros(N,H-1); per_rate_body = zeros(N,H-1);

z = model.R'*yc;              % latent in standardized space
for n = 1:N
    k0 = starts(n);
    zk = z(:,k0);
    ypred = zeros(p,H);
    for h = 1:H
        zk = model.A*zk + model.B*uc(:,k0+h-1);
        ypred(:,h) = model.y_mean + model.P*zk;
    end
    ypred_raw = ypred .* scales.y_scale + scales.y_offset;
    ytrue = y_raw(:,k0+1:k0+H);
    per_raw = repmat(y_raw(:,k0),1,H);   % persistence = hold last observation

    for h = 1:H
        err_pos(n,h) = norm(ypred_raw(1:3,h)-ytrue(1:3,h));
        err_att(n,h) = norm(ypred_raw(4:6,h)-ytrue(4:6,h));
        per_pos(n,h) = norm(per_raw(1:3,h)-ytrue(1:3,h));
        per_att(n,h) = norm(per_raw(4:6,h)-ytrue(4:6,h));
    end
    % Speed: dataset convention v(k) = 100*(pos(k+1)-pos(k)) in m/s.
    vpred = 100*diff(ypred_raw(1:3,:),1,2);
    vtrue = 100*diff(ytrue(1:3,:),1,2);
    erate_pred = 100*diff(ypred_raw(4:6,:),1,2);
    erate_true = 100*diff(ytrue(4:6,:),1,2);
    for h = 1:H-1
        err_vel(n,h)  = norm(vpred(:,h)-vtrue(:,h));
        per_vel(n,h)  = norm(vtrue(:,h)); % persistence speed = 0
        err_rate(n,h) = norm(erate_pred(:,h)-erate_true(:,h));
        per_rate(n,h) = norm(erate_true(:,h)); % persistence rate = 0
        % body-frame rate: pqr = T(att_true) * euler_rate, standard ZYX
        phi=ytrue(1,h+1); th=ytrue(2,h+1); % use true attitude for the transform
        Tatt=[1 0 -sin(th); 0 cos(phi) sin(phi)*cos(th); 0 -sin(phi) cos(phi)*cos(th)];
        err_rate_body(n,h) = norm(Tatt*(erate_pred(:,h)-erate_true(:,h)));
        per_rate_body(n,h) = norm(Tatt*erate_true(:,h));
    end
end

out=struct('err_pos',err_pos,'err_att',err_att,'err_vel',err_vel, ...
    'err_rate',err_rate,'err_rate_body',err_rate_body, ...
    'per_pos',per_pos,'per_att',per_att,'per_vel',per_vel, ...
    'per_rate',per_rate,'per_rate_body',per_rate_body, ...
    'starts',starts,'N',N,'H',H,'fs',100,'horizon_s',H/100);
end
