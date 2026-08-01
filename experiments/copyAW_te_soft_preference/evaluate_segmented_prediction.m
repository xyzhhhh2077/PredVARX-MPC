function out = evaluate_segmented_prediction(y,u,run_id,Etask,model)
%EVALUATE_SEGMENTED_PREDICTION One-step validation without crossing run boundaries.
run_id=run_id(:)'; valid=find(run_id(1:end-1)==run_id(2:end));
yc=y-model.y_mean; uc=u-model.u_mean;
z=model.R'*yc;
z_actual=z(:,valid+1);
z_predicted=model.A*z(:,valid)+model.B*uc(:,valid);
z_persistence=z(:,valid);
y_predicted=model.y_mean+model.P*z_predicted;
task_actual=Etask'*y(:,valid+1);
task_predicted=Etask'*y_predicted;
task_persistence=Etask'*y(:,valid);
task_error=task_actual-task_predicted;
den=sum((task_actual-mean(task_actual,2)).^2,2);
out=struct();
out.transition_index=valid;
out.run_id=run_id(valid+1);
out.transition_count=numel(valid);
out.latent_actual=z_actual;
out.latent_predicted=z_predicted;
out.latent_rmse=sqrt(mean((z_actual-z_predicted).^2,2));
out.latent_persistence_rmse=sqrt(mean((z_actual-z_persistence).^2,2));
out.task_actual=task_actual;
out.task_predicted=task_predicted;
out.task_rmse=sqrt(mean(task_error.^2,2));
out.task_persistence_rmse=sqrt(mean((task_actual-task_persistence).^2,2));
out.task_r2=1-sum(task_error.^2,2)./max(den,eps);
out.full_output_rmse=sqrt(mean((y(:,valid+1)-y_predicted).^2,'all'));
end