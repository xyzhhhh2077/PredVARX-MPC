%% copyAX CDR diagnostic verification (model -> parameters -> results)
cd('C:/Users/ROG/predvarx-repo/experiments/copyAX_cdr_soft_preference');
addpath(pwd);
D = load(fullfile(pwd,'data','copyAX_cdr_dataset.mat'));
S = load(fullfile(pwd,'results','copyAX_cdr_soft_preference_data.mat'));

segment_id = double(D.segment_id(:)');
train_mask = ismember(segment_id,1:6);
validation_mask = ismember(segment_id,7:8);
ytr = double(D.y(:,train_mask)); utr = double(D.u(:,train_mask));
yva = double(D.y(:,validation_mask)); uva = double(D.u(:,validation_mask));
trseg = segment_id(train_mask); vaseg = segment_id(validation_mask);
yoff = mean(ytr,2); ys = std(ytr,0,2); uoff = mean(utr,2); us = std(utr,0,2);
ytr = (ytr-yoff)./ys; utr = (utr-uoff)./us; yva = (yva-yoff)./ys; uva = (uva-uoff)./us;

fprintf('== MODEL ADEQUACY ==\n');
fprintf('spectral_radius %.6f (true DC mode exp(+0.001) per step: mildly growing)\n', ...
    S.fitstats.spectral_radius);
fprintf('dual_error %.3e   anchor_preservation %.3e\n', ...
    S.fitstats.dual_error, S.fitstats.anchor_preservation_error);

valid = find(trseg(1:end-1)==trseg(2:end));
z = S.model.R'*ytr;
zn = z(:,valid+1); zc = z(:,valid); ur = utr(:,valid);
Eps = zn - S.model.A*zc - S.model.B*ur;
fprintf('residual std by latent dim: %s\n', mat2str(std(Eps,0,2)',4));
fprintf('residual autocorrelation (max |acf| over 5 dims):\n');
for lag = 1:5
    ac = zeros(5,1);
    for d = 1:5
        e1 = Eps(d,lag+1:end); e2 = Eps(d,1:end-lag);
        ac(d) = abs(sum(e1.*e2))/sqrt(sum(e1.^2)*sum(e2.^2));
    end
    fprintf('  lag %d: %.3f\n', lag, max(ac));
end

fprintf('\n== PARAMETER ADEQUACY ==\n');
scale = max(trace(z*z'/size(z,2))/5,1e-12);
G = PhiGram(zc,ur,1e-8*max(scale,1));
fprintf('regression Gram cond: %.3e (65 params over %d transitions, ratio %.1f)\n', ...
    cond(G), numel(valid), numel(valid)/65);
ev = sort(real(eig(S.fitstats.S_yu)),'descend');
fprintf('S_yu (input-driven latent cov) top eigenvalues: %s\n', mat2str(ev(1:min(8,numel(ev)))',4));
fprintf('CRTE eigenvalues (free block): %s\n', mat2str(S.fitstats.crte_eigenvalues',4));
fprintf('cond(C_y): %.3e\n', cond(S.fitstats.C_y));

fprintf('\n== RESULT VALIDITY ==\n');
validv = find(vaseg(1:end-1)==vaseg(2:end));
zva = S.model.R'*yva;
zact = zva(:,validv+1); zpre = S.model.A*zva(:,validv)+S.model.B*uva(:,validv);
ypre = S.model.y_mean + S.model.P*zpre;
taskact = S.Etask'*yva(:,validv+1); taskpre = S.Etask'*ypre;
for seg = [7 8]
    selm = vaseg(validv+1)==seg;          % logical mask over transitions j=1..2000
    err = taskact(:,selm)-taskpre(:,selm);
    den = sum((taskact(:,selm)-mean(taskact(:,selm),2)).^2,2);
    r2 = 1 - sum(err.^2,2)./max(den,eps);
    fprintf('validation segment %d: task R2 = %s (n=%d)\n', seg, mat2str(r2',4), nnz(selm));
end
fprintf('contribution top sensors: %s\n', mat2str(topk(S.dirstats.contribution,6)',4));
fprintf('contribution min/max: %.4f / %.4f\n', ...
    min(S.dirstats.contribution), max(S.dirstats.contribution));

function G = PhiGram(zc,ur,ridge)
G = [zc;ur]*[zc;ur]' + ridge*eye(size(zc,1)+size(ur,1));
G = (G+G')/2;
end

function v = topk(x,k)
[~,idx] = sort(x,'descend'); v = idx(1:k);
end
