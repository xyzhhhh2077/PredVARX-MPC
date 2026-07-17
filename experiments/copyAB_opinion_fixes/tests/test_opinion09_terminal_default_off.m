function test_opinion09_terminal_default_off
% Opinion 9: optional terminal cost defaults OFF and preserves old call shape.
% Strict convexity of the QP is not recursive feasibility and not stability.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

A = 0.8; B = 1; P = 1; R = 1;
model.A = A; model.B = B; model.P = P; model.R = R;
model.y_mean = 2; model.u_mean = 3;
model.Sigma_eps = 0.04; model.Sigma_obs = 0.01;

opt.N = 3; opt.Q = 10; opt.Ru = 0.1;
opt.u_min = -10; opt.u_max = 10;
opt.H = [1; -1]; opt.h = [5; 5]; opt.alpha_joint = 0.10;
% Intentionally omit use_terminal_cost → default OFF.

[z, y_pred, U, out] = centered_smpc_step(2, 2, model, opt);

assert(isnumeric(z) && isscalar(z), 'z must be numeric scalar for 1-D model');
assert(isnumeric(y_pred) && ~isempty(y_pred), 'y_pred must be non-empty numeric');
assert(isnumeric(U) && isequal(size(U), [opt.N*size(model.B,2), 1]), ...
    'U must be N*nu x 1');
assert(isstruct(out), 'out must be a struct');
assert(isfield(out,'cost') && isfield(out,'exitflag') && isfield(out,'A_ch'), ...
    'out must keep legacy diagnostic fields');
assert(isfield(out,'use_terminal_cost') && out.use_terminal_cost == false, ...
    'default must leave use_terminal_cost false');
assert(isfield(out,'terminal_cost_applied') && out.terminal_cost_applied == false, ...
    'default must not apply terminal cost');
assert(isempty(out.Pterm), 'default Pterm must be empty');
assert(out.exitflag > 0, 'default QP must be feasible');
assert(all(U >= opt.u_min-1e-10 & U <= opt.u_max+1e-10), 'input bounds');
assert(max(out.A_ch*U - out.b_ch) < 1e-7, 'chance rows must hold under default');

% Explicit false must match omitted-field path (same-shaped plan).
opt_off = opt;
opt_off.use_terminal_cost = false;
[~, ~, U_off, out_off] = centered_smpc_step(2, 2, model, opt_off);
assert(norm(U - U_off) < 1e-10, 'explicit false must match default-off plan');
assert(abs(out.cost - out_off.cost) < 1e-10, 'explicit false must match default cost');
assert(out_off.terminal_cost_applied == false);

% Optional ON path is callable; if dlyap succeeds it should mark applied.
opt_on = opt;
opt_on.use_terminal_cost = true;
try
    [~, ~, U_on, out_on] = centered_smpc_step(2, 2, model, opt_on);
    assert(isfield(out_on,'terminal_cost_applied'));
    if out_on.terminal_cost_applied
        assert(~isempty(out_on.Pterm) && all(isfinite(out_on.Pterm(:))), ...
            'applied Pterm must be finite');
        % Terminal weight should generally change the plan vs default-off.
        % (not a stability claim — only a non-degeneracy check)
        assert(norm(U_on - U) > 1e-12 || abs(out_on.cost - out.cost) > 1e-12, ...
            'terminal cost applied but objective/plan identical');
    end
catch ME
    error('use_terminal_cost=true path must remain callable: %s', ME.message);
end

fprintf(['PASS opinion09: default terminal OFF preserves call shape; ' ...
    'optional ON callable (applied=%d).\n'], logical(out_on.terminal_cost_applied));
end
