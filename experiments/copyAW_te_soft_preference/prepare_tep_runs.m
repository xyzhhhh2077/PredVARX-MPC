function data = prepare_tep_runs(raw,train_runs,validation_runs)
%PREPARE_TEP_RUNS Split the Rieth TEP table by complete simulation runs.
train_runs = unique(train_runs(:)');
validation_runs = unique(validation_runs(:)');
if any(ismember(train_runs,validation_runs))
    error('prepare_tep_runs:OverlappingRuns', ...
        'Training and validation simulation runs must be disjoint.');
end

required = [{'faultNumber','simulationRun','sample'}, ...
    arrayfun(@(j)sprintf('xmeas_%d',j),1:41,'UniformOutput',false), ...
    arrayfun(@(j)sprintf('xmv_%d',j),1:11,'UniformOutput',false)];
missing = setdiff(required,raw.Properties.VariableNames);
assert(isempty(missing),'TEP table is missing required columns: %s',strjoin(missing,', '));
assert(all(raw.faultNumber==0),'Only fault-free TEP runs are accepted for this experiment.');

train = select_runs(raw,train_runs);
validation = select_runs(raw,validation_runs);
measurement_names = arrayfun(@(j)sprintf('xmeas_%d',j),1:41,'UniformOutput',false);
input_names = arrayfun(@(j)sprintf('xmv_%d',j),1:11,'UniformOutput',false);

data = struct();
data.y_train = table2array(train(:,measurement_names))';
data.u_train = table2array(train(:,input_names))';
data.train_run_id = train.simulationRun';
data.train_sample = train.sample';
data.y_validation = table2array(validation(:,measurement_names))';
data.u_validation = table2array(validation(:,input_names))';
data.validation_run_id = validation.simulationRun';
data.validation_sample = validation.sample';
data.measurement_names = measurement_names;
data.input_names = input_names;
end

function selected = select_runs(raw,runs)
selected = raw(ismember(raw.simulationRun,runs),:);
assert(~isempty(selected),'Requested TEP run split is empty.');
selected = sortrows(selected,{'simulationRun','sample'});
found = unique(selected.simulationRun)';
assert(isequal(found,runs),'One or more requested TEP runs are absent.');
end