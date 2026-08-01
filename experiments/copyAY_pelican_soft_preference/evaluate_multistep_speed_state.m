function out = evaluate_multistep_speed_state(y_raw,u_raw,run_id,model,scales,H,start_step)
%EVALUATE_MULTISTEP_SPEED_STATE Recursive multi-step prediction on speed state.
% y_raw = [Vel(3); pqr(3); Motors(4)] in physical units (m/s, rad/s, dimless).
% The ICRA 2018 benchmark predicts translational velocity and body rates
% directly -- these ARE the output channels, so per-horizon errors are
% measured without any integration/differencing:
%   err_vel(n,h) = ||Vel_pred(:,h) - Vel_true(:,h)||       (m/s)
%   err_rate(n,h) = ||pqr_pred(:,h) - pqr_true(:,h)||      (rad/s)
% Persistence baseline holds the last observation (vel/rate stay constant).
% Outputs per start point (n) and horizon (h):
%   err_vel, err_rate, per_vel, per_rate   N x (H-1)
%   starts  1 x N  start indices used
%   N        number of starts
if nargin < 7 || isempty(start_step), start_step = H; end
[p,T] = size(y_raw); m = size(u_raw,1);
run_id = run_id(:)';
assert(size(u_raw,2)==T && numel(run_id)==T,'Sample counts must match.');
assert(p==10,'Expected Pelican speed state (10 channels).');
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

err_vel = zeros(N,H-1); err_rate = zeros(N,H-1);
per_vel = zeros(N,H-1); per_rate = zeros(N,H-1);
mae_vel = zeros(N,H-1); mae_rate = zeros(N,H-1);
mae_pvel = zeros(N,H-1); mae_prate = zeros(N,H-1);

z = model.R'*yc;              % latent in standardized space
for n = 1:N
    k0 = starts(n);
    zk = z(:,k0);
    ypred = zeros(p,H-1);
    % The first predicted step uses z_{k0+1} = A z_{k0} + B u_{k0}, i.e.
    % y_true(:,k0+1) is the first target (standard causality contract).
    for h = 1:H-1
        zk = model.A*zk + model.B*uc(:,k0+h-1);
        ypred(:,h) = model.y_mean + model.P*zk;
    end
    ypred_raw = ypred .* scales.y_scale + scales.y_offset;
    ytrue = y_raw(:,k0+1:k0+H-1);
    per_raw = repmat(y_raw(:,k0),1,H-1);   % persistence = hold last observation

    err_vel(n,:)  = sqrt(sum((ypred_raw(1:3,:)-ytrue(1:3,:)).^2,1));
    err_rate(n,:) = sqrt(sum((ypred_raw(4:6,:)-ytrue(4:6,:)).^2,1));
    per_vel(n,:)  = sqrt(sum((per_raw(1:3,:)-ytrue(1:3,:)).^2,1));
    per_rate(n,:) = sqrt(sum((per_raw(4:6,:)-ytrue(4:6,:)).^2,1));
    % paper norm (30): mean absolute error per axis
    mae_vel(n,:)  = sum(abs(ypred_raw(1:3,:)-ytrue(1:3,:)),1)/3;
    mae_rate(n,:) = sum(abs(ypred_raw(4:6,:)-ytrue(4:6,:)),1)/3;
    mae_pvel(n,:) = sum(abs(per_raw(1:3,:)-ytrue(1:3,:)),1)/3;
    mae_prate(n,:)= sum(abs(per_raw(4:6,:)-ytrue(4:6,:)),1)/3;
end

out=struct('err_vel',err_vel,'err_rate',err_rate, ...
    'per_vel',per_vel,'per_rate',per_rate, ...
    'mae_vel',mae_vel,'mae_rate',mae_rate, ...
    'mae_pvel',mae_pvel,'mae_prate',mae_prate, ...
    'starts',starts,'N',N,'H',H,'fs',100,'horizon_s',H/100);
end
