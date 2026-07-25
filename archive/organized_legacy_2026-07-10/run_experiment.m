function run_experiment(script_path)
%RUN_EXPERIMENT Run a copied MATLAB experiment in its own directory.
%   It preserves the caller's working folder and makes the canonical
%   PredVARX identifier available without modifying historical scripts.

    root = fileparts(mfilename('fullpath'));
    id_dir = fullfile(root, '02_identification');
    original_dir = pwd;
    cleanup = onCleanup(@() cd(original_dir)); %#ok<NASGU>

    if nargin ~= 1 || ~(ischar(script_path) || isstring(script_path))
        error('Usage: run_experiment(full_path_to_script)');
    end

    script_path = char(script_path);
    if ~isfile(script_path)
        error('Script not found: %s', script_path);
    end

    addpath(id_dir, '-begin');
    cd(fileparts(script_path));
    run(script_path);
end
