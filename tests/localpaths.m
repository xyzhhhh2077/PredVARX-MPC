function p = localpaths()
%LOCALPATHS Canonical source roots for PredVARX-MPC unit tests.
here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
p.repo = repo;
p.main = fullfile(repo, 'main');
p.copyO = fullfile(repo, 'experiments', 'copyO_oblique');
p.copyP = fullfile(repo, 'experiments', 'copyP_centered_smpc');
p.copyQ = fullfile(repo, 'experiments', 'copyQ_control_aware');
p.copyR = fullfile(repo, 'experiments', 'copyR_moqin_oblique');
p.copyAA = fullfile(repo, 'experiments', 'copyAA_split_control_free_oblique');
end
