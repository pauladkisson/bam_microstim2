function mse_loss = eval_thresh_cor(I_els, cable_thresholds, thresh_cor, depols)
    z_bounds = [0, 2000]*1e-4;
    z_res = 1e-5;
    lif_thresholds = zeros(length(I_els), length(depols));
    mse_loss = 0;
    for i = 1:length(I_els)
        I_el = I_els(i);
        for j = 1:length(depols)
            depol = depols(j);
            mirror_est = @(z) mir_est(z, I_el, thresh_cor, depol);
            z_thresh = bisect_search(mirror_est, @lif_eval, z_bounds, z_res);
            lif_thresholds(i, j) = z_thresh;
            loss = (z_thresh - cable_thresholds(i, j))^2;
            mse_loss = mse_loss + loss;
        end
    end
    mse_loss = mse_loss / (length(I_els)*length(depols));
end