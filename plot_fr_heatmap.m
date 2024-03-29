%%% Paul Adkisson
%%% 2/14/2023
%%% Plot Mean Firing Rate Heatmaps (over space) for each Stimulation Condition
function plot_fr_heatmap(sim_name, pulse_amps, stim_amps, t, t_cut, ex_c, ...
    pulse_coherences, galvanic_coherences, control_coherences, anodic_coherences, ...
    default_colors, start_trial, end_trial, num_trials, num_group, num_affected, plot_name)
    if plot_name == "p1_wins"
        win_num = 1;
    else
        win_num = 2;
    end
    dt = t(2) - t(1);
    win_size = 50e-3;
    avg_win_size = 200e-3;
    n_pts_per_bin = 4; %Number of time points per averaging window
    down_factor = floor(avg_win_size/(n_pts_per_bin*dt));
    t_down = t(1:down_factor:end);
    stim_frs = zeros(length(stim_amps), num_trials, length(t_down), num_group);
    stim_cvs = zeros(length(stim_amps), num_trials, length(t_down));
    stim_skews = zeros(length(stim_amps), num_trials, length(t_down));
    stim_kurts = zeros(length(stim_amps), num_trials, length(t_down));
    load(sprintf("Simulation %s/ustim/r.mat", sim_name), "ball_r");
    ball_rs = get_ball_rs(ball_r, num_affected, num_group);
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        c = ex_c(j);
        pulse = j<=length(pulse_amps);
        if pulse
            disp("Pulsatile")
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
        [sim_name, stim_amp*1e6]);
        stim_coherences = pulse_coherences;
        else
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
        [sim_name, stim_amp*1e6]);
            if stim_amp < 0 %cathodic GS
                disp("Cathodic GS")
                stim_coherences = galvanic_coherences;
            elseif stim_amp == 0 % control
                disp("Control")
                stim_coherences = control_coherences;
            else
                disp("Anodic GS")
                stim_coherences = anodic_coherences;
            end
        end
        load(strcat(output_stimpath, "/decisions.mat"), "decisions", "decision_times")
        for trial = start_trial:end_trial
            fprintf("Trial: %0.0f \n", trial)
            relative_trial = trial - start_trial + 1;
            dec_time = decision_times(relative_trial, stim_coherences==c);
            if decisions(relative_trial, stim_coherences==c) ~= win_num || ...
                    dec_time > t_cut
                stim_frs(j, relative_trial, :) = NaN;
                stim_cvs(j, relative_trial, :) = NaN;
                stim_skews(j, relative_trial, :) = NaN;
                stim_kurts(j, relative_trial, :) = NaN;
                continue %skip trials where P1 doesn't win/lose or decision takes too long
            end
            load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                   "recspikes")
            neuron_frs = recspikes2neuron_frs(recspikes, win_size, avg_win_size, t, num_group);
            % Downsample for appropriate time-frequency resolution
            neuron_frs = neuron_frs(1:down_factor:end, :);
            stim_frs(j, trial, :, :) = neuron_frs;
            stim_cvs(j, trial, :) = std(neuron_frs, [], 2) ./ mean(neuron_frs, 2);
            stim_skews(j, trial, :) = skewness(neuron_frs, 0, 2); %flag = 0 to correct for sample bias
            stim_kurts(j, trial, :) = kurtosis(neuron_frs, 0, 2); %flag = 0 to correct for sample bias
        end
    end
    ps_frs = reshape(mean(stim_frs(1, :, :, :), 2, 'omitnan'), [length(t_down), num_group]);
    ps_cv = reshape(mean(stim_cvs(1, :, :), 2, 'omitnan'), [length(t_down), 1]);
    ps_skew = reshape(mean(stim_skews(1, :, :), 2, 'omitnan'), [length(t_down), 1]);
    ps_kurt = reshape(mean(stim_kurts(1, :, :), 2, 'omitnan'), [length(t_down), 1]);
    
    gs_frs = reshape(mean(stim_frs(2, :, :, :), 2, 'omitnan'), [length(t_down), num_group]);
    gs_cv = reshape(mean(stim_cvs(2, :, :), 2, 'omitnan'), [length(t_down), 1]);
    gs_skew = reshape(mean(stim_skews(2, :, :), 2, 'omitnan'), [length(t_down), 1]);
    gs_kurt = reshape(mean(stim_kurts(2, :, :), 2, 'omitnan'), [length(t_down), 1]);
    
    ctrl_frs = reshape(mean(stim_frs(3, :, :, :), 2, 'omitnan'), [length(t_down), num_group]);
    ctrl_cv = reshape(mean(stim_cvs(3, :, :), 2, 'omitnan'), [length(t_down), 1]);
    ctrl_skew = reshape(mean(stim_skews(3, :, :), 2, 'omitnan'), [length(t_down), 1]);
    ctrl_kurt = reshape(mean(stim_kurts(3, :, :), 2, 'omitnan'), [length(t_down), 1]);
    
    an_frs = reshape(mean(stim_frs(4, :, :, :), 2, 'omitnan'), [length(t_down), num_group]);
    an_cv = reshape(mean(stim_cvs(4, :, :), 2, 'omitnan'), [length(t_down), 1]);
    an_skew = reshape(mean(stim_skews(4, :, :), 2, 'omitnan'), [length(t_down), 1]);
    an_kurt = reshape(mean(stim_kurts(4, :, :), 2, 'omitnan'), [length(t_down), 1]);
    
    % Statistics
    stim_trials = [sum(~isnan(stim_frs(1, :, 1, 1))), ...
                   sum(~isnan(stim_frs(2, :, 1, 1))), ...
                   sum(~isnan(stim_frs(3, :, 1, 1))), ...
                   sum(~isnan(stim_frs(4, :, 1, 1)))  ];
    stim_skew_ses = skew_se(stim_trials);
    stim_kurt_ses = kurt_se(stim_trials);
    beginning_of_task = t_down>=1.1 & t_down<=1.2;
    end_of_task = t_down>=2.9 & t_down<=3;
    
    ps_cv_start = reshape(mean(stim_cvs(1, :, beginning_of_task), 3), [num_trials, 1]);
    gs_cv_start = reshape(mean(stim_cvs(2, :, beginning_of_task), 3), [num_trials, 1]);
    ctrl_cv_start = reshape(mean(stim_cvs(3, :, beginning_of_task), 3), [num_trials, 1]);
    an_cv_start = reshape(mean(stim_cvs(4, :, beginning_of_task), 3), [num_trials, 1]);
    ps_cv_end = reshape(mean(stim_cvs(1, :, end_of_task), 3), [num_trials, 1]);
    gs_cv_end = reshape(mean(stim_cvs(2, :, end_of_task), 3), [num_trials, 1]);
    ctrl_cv_end = reshape(mean(stim_cvs(3, :, end_of_task), 3), [num_trials, 1]);
    an_cv_end = reshape(mean(stim_cvs(4, :, end_of_task), 3), [num_trials, 1]);
    
    ps_skew_start = reshape(mean(stim_skews(1, :, beginning_of_task), 3), [num_trials, 1]);
    gs_skew_start = reshape(mean(stim_skews(2, :, beginning_of_task), 3), [num_trials, 1]);
    ctrl_skew_start = reshape(mean(stim_skews(3, :, beginning_of_task), 3), [num_trials, 1]);
    an_skew_start = reshape(mean(stim_skews(4, :, beginning_of_task), 3), [num_trials, 1]);
    ps_skew_end = reshape(mean(stim_skews(1, :, end_of_task), 3), [num_trials, 1]);
    gs_skew_end = reshape(mean(stim_skews(2, :, end_of_task), 3), [num_trials, 1]);
    ctrl_skew_end = reshape(mean(stim_skews(3, :, end_of_task), 3), [num_trials, 1]);
    an_skew_end = reshape(mean(stim_skews(4, :, end_of_task), 3), [num_trials, 1]);
    
    ps_kurt_start = reshape(mean(stim_kurts(1, :, beginning_of_task), 3), [num_trials, 1]);
    gs_kurt_start = reshape(mean(stim_kurts(2, :, beginning_of_task), 3), [num_trials, 1]);
    ctrl_kurt_start = reshape(mean(stim_kurts(3, :, beginning_of_task), 3), [num_trials, 1]);
    an_kurt_start = reshape(mean(stim_kurts(4, :, beginning_of_task), 3), [num_trials, 1]);
    ps_kurt_end = reshape(mean(stim_kurts(1, :, end_of_task), 3), [num_trials, 1]);
    gs_kurt_end = reshape(mean(stim_kurts(2, :, end_of_task), 3), [num_trials, 1]);
    ctrl_kurt_end = reshape(mean(stim_kurts(3, :, end_of_task), 3), [num_trials, 1]);
    an_kurt_end = reshape(mean(stim_kurts(4, :, end_of_task), 3), [num_trials, 1]);
    
%     ps_skew_start = mean(ps_skew(beginning_of_task));
%     gs_skew_start = mean(gs_skew(beginning_of_task));
%     ctrl_skew_start = mean(ctrl_skew(beginning_of_task)); 
%     an_skew_start = mean(an_skew(beginning_of_task));
%     ps_skew_end = mean(ps_skew(end_of_task));
%     gs_skew_end = mean(gs_skew(end_of_task));
%     ctrl_skew_end = mean(ctrl_skew(end_of_task)); 
%     an_skew_end = mean(an_skew(end_of_task));
%     
%     ps_kurt_start = mean(ps_kurt(beginning_of_task));
%     gs_kurt_start = mean(gs_kurt(beginning_of_task));
%     ctrl_kurt_start = mean(ctrl_kurt(beginning_of_task)); 
%     an_kurt_start = mean(an_kurt(beginning_of_task));
%     ps_kurt_end = mean(ps_kurt(end_of_task));
%     gs_kurt_end = mean(gs_kurt(end_of_task));
%     ctrl_kurt_end = mean(ctrl_kurt(end_of_task)); 
%     an_kurt_end = mean(an_kurt(end_of_task));
    
    [~, ~, stats] = kruskalwallis([ps_cv_start, gs_cv_start, ctrl_cv_start, an_cv_start], [], 'off');
    c = multcompare(stats, 'Display', 'off');
    p_ps_cgs_cv_start = c(1, end);
    p_ags_ctrl_cv_start = c(end, end);
    ps_cv_start_quantiles = quantile(ps_cv_start, [0.25, 0.5, 0.75]);
    cgs_cv_start_quantiles = quantile(gs_cv_start, [0.25, 0.5, 0.75]);
    ags_cv_start_quantiles = quantile(an_cv_start, [0.25, 0.5, 0.75]);
    ctrl_cv_start_quantiles = quantile(ctrl_cv_start, [0.25, 0.5, 0.75]);
    
    [~, ~, stats] = kruskalwallis([ps_cv_end, gs_cv_end, ctrl_cv_end, an_cv_end], [], 'off');
    c = multcompare(stats, 'Display', 'off');
    p_ps_cgs_cv_end = c(1, end);
    p_ags_ctrl_cv_end = c(end, end);
    ps_cv_end_quantiles = quantile(ps_cv_end, [0.25, 0.5, 0.75]);
    cgs_cv_end_quantiles = quantile(gs_cv_end, [0.25, 0.5, 0.75]);
    ags_cv_end_quantiles = quantile(an_cv_end, [0.25, 0.5, 0.75]);
    ctrl_cv_end_quantiles = quantile(ctrl_cv_end, [0.25, 0.5, 0.75]);
    
    [~, ~, stats] = kruskalwallis([ps_skew_start, gs_skew_start, ctrl_skew_start, an_skew_start], [], 'off');
    c = multcompare(stats, 'Display', 'off');
    p_ps_cgs_skew_start = c(1, end);
    p_ags_ctrl_skew_start = c(end, end);
    ps_skew_start_quantiles = quantile(ps_skew_start, [0.25, 0.5, 0.75]);
    cgs_skew_start_quantiles = quantile(gs_skew_start, [0.25, 0.5, 0.75]);
    ags_skew_start_quantiles = quantile(an_skew_start, [0.25, 0.5, 0.75]);
    ctrl_skew_start_quantiles = quantile(ctrl_skew_start, [0.25, 0.5, 0.75]);
    
    [~, ~, stats] = kruskalwallis([ps_skew_end, gs_skew_end, ctrl_skew_end, an_skew_end], [], 'off');
    c = multcompare(stats, 'Display', 'off');
    p_ps_cgs_skew_end = c(1, end);
    p_ags_ctrl_skew_end = c(end, end);
    ps_skew_end_quantiles = quantile(ps_skew_end, [0.25, 0.5, 0.75]);
    cgs_skew_end_quantiles = quantile(gs_skew_end, [0.25, 0.5, 0.75]);
    ags_skew_end_quantiles = quantile(an_skew_end, [0.25, 0.5, 0.75]);
    ctrl_skew_end_quantiles = quantile(ctrl_skew_end, [0.25, 0.5, 0.75]);
    
    [~, ~, stats] = kruskalwallis([ps_kurt_start, gs_kurt_start, ctrl_kurt_start, an_kurt_start], [], 'off');
    c = multcompare(stats, 'Display', 'off');
    p_ps_cgs_kurt_start = c(1, end);
    p_ags_ctrl_kurt_start = c(end, end);
    ps_kurt_start_quantiles = quantile(ps_kurt_start, [0.25, 0.5, 0.75]);
    cgs_kurt_start_quantiles = quantile(gs_kurt_start, [0.25, 0.5, 0.75]);
    ags_kurt_start_quantiles = quantile(an_kurt_start, [0.25, 0.5, 0.75]);
    ctrl_kurt_start_quantiles = quantile(ctrl_kurt_start, [0.25, 0.5, 0.75]);
    
    [~, ~, stats] = kruskalwallis([ps_kurt_end, gs_kurt_end, ctrl_kurt_end, an_kurt_end], [], 'off');
    c = multcompare(stats, 'Display', 'off');
    p_ps_cgs_kurt_end = c(1, end);
    p_ags_ctrl_kurt_end = c(end, end);
    ps_kurt_end_quantiles = quantile(ps_kurt_end, [0.25, 0.5, 0.75]);
    cgs_kurt_end_quantiles = quantile(gs_kurt_end, [0.25, 0.5, 0.75]);
    ags_kurt_end_quantiles = quantile(an_kurt_end, [0.25, 0.5, 0.75]);
    ctrl_kurt_end_quantiles = quantile(ctrl_kurt_end, [0.25, 0.5, 0.75]);
    
    % Start of Trial
    fprintf([...
            'start-of-trial: PS induced a different CV (%0.2f, %0.2f, %0.2f)', ...
            ' than CGS (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ps_cv_start_quantiles(1), ps_cv_start_quantiles(2), ps_cv_start_quantiles(3), ...
            cgs_cv_start_quantiles(1), cgs_cv_start_quantiles(2), cgs_cv_start_quantiles(3), p_ps_cgs_cv_start)
    fprintf([...
            'start-of-trial: AGS induced a equivalent CV (%0.2f, %0.2f, %0.2f)', ...
            ' as Control (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ags_cv_start_quantiles(1), ags_cv_start_quantiles(2), ags_cv_start_quantiles(3), ...
            ctrl_cv_start_quantiles(1), ctrl_cv_start_quantiles(2), ctrl_cv_start_quantiles(3), p_ags_ctrl_cv_start)
        
    fprintf([...
            'start-of-trial: PS induced a different Skew (%0.2f, %0.2f, %0.2f)', ...
            ' than CGS (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ps_skew_start_quantiles(1), ps_skew_start_quantiles(2), ps_skew_start_quantiles(3), ...
            cgs_skew_start_quantiles(1), cgs_skew_start_quantiles(2), cgs_skew_start_quantiles(3), p_ps_cgs_skew_start)
    fprintf([...
            'start-of-trial: AGS induced a equivalent Skew (%0.2f, %0.2f, %0.2f)', ...
            ' as Control (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ags_skew_start_quantiles(1), ags_skew_start_quantiles(2), ags_skew_start_quantiles(3), ...
            ctrl_skew_start_quantiles(1), ctrl_skew_start_quantiles(2), ctrl_skew_start_quantiles(3), p_ags_ctrl_skew_start)
        
    fprintf([...
            'start-of-trial: PS induced a different Kurtosis (%0.2f, %0.2f, %0.2f)', ...
            ' than CGS (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ps_kurt_start_quantiles(1), ps_kurt_start_quantiles(2), ps_kurt_start_quantiles(3), ...
            cgs_kurt_start_quantiles(1), cgs_kurt_start_quantiles(2), cgs_kurt_start_quantiles(3), p_ps_cgs_kurt_start)
    fprintf([...
            'start-of-trial: AGS induced a equivalent Kurtosis (%0.2f, %0.2f, %0.2f)', ...
            ' as Control (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ags_kurt_start_quantiles(1), ags_kurt_start_quantiles(2), ags_kurt_start_quantiles(3), ...
            ctrl_kurt_start_quantiles(1), ctrl_kurt_start_quantiles(2), ctrl_kurt_start_quantiles(3), p_ags_ctrl_kurt_start)
        
    
    % End of Trial
    fprintf([...
            'end-of-trial: PS induced a different CV (%0.2f, %0.2f, %0.2f)', ...
            ' than CGS (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ps_cv_end_quantiles(1), ps_cv_end_quantiles(2), ps_cv_end_quantiles(3), ...
            cgs_cv_end_quantiles(1), cgs_cv_end_quantiles(2), cgs_cv_end_quantiles(3), p_ps_cgs_cv_end)
    fprintf([...
            'end-of-trial: AGS induced a equivalent CV (%0.2f, %0.2f, %0.2f)', ...
            ' as Control (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ags_cv_end_quantiles(1), ags_cv_end_quantiles(2), ags_cv_end_quantiles(3), ...
            ctrl_cv_end_quantiles(1), ctrl_cv_end_quantiles(2), ctrl_cv_end_quantiles(3), p_ags_ctrl_cv_end)
        
    fprintf([...
            'end-of-trial: PS induced a different Skew (%0.2f, %0.2f, %0.2f)', ...
            ' than CGS (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ps_skew_end_quantiles(1), ps_skew_end_quantiles(2), ps_skew_end_quantiles(3), ...
            cgs_skew_end_quantiles(1), cgs_skew_end_quantiles(2), cgs_skew_end_quantiles(3), p_ps_cgs_skew_end)
    fprintf([...
            'end-of-trial: AGS induced a equivalent Skew (%0.2f, %0.2f, %0.2f)', ...
            ' as Control (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ags_skew_end_quantiles(1), ags_skew_end_quantiles(2), ags_skew_end_quantiles(3), ...
            ctrl_skew_end_quantiles(1), ctrl_skew_end_quantiles(2), ctrl_skew_end_quantiles(3), p_ags_ctrl_skew_end)
        
    fprintf([...
            'end-of-trial: PS induced a different Kurtosis (%0.2f, %0.2f, %0.2f)', ...
            ' than CGS (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ps_kurt_end_quantiles(1), ps_kurt_end_quantiles(2), ps_kurt_end_quantiles(3), ...
            cgs_kurt_end_quantiles(1), cgs_kurt_end_quantiles(2), cgs_kurt_end_quantiles(3), p_ps_cgs_kurt_end)
    fprintf([...
            'end-of-trial: AGS induced a equivalent Kurtosis (%0.2f, %0.2f, %0.2f)', ...
            ' as Control (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'],  ...
            ags_kurt_end_quantiles(1), ags_kurt_end_quantiles(2), ags_kurt_end_quantiles(3), ...
            ctrl_kurt_end_quantiles(1), ctrl_kurt_end_quantiles(2), ctrl_kurt_end_quantiles(3), p_ags_ctrl_kurt_end)
   
%     fprintf([...
%         'start-of-trial: PS   induced increases in cv (%0.2f +/- %0.2f), '...
%         'skew (%0.2f +/- %0.2f) and kurtosis (%0.2f +/- %0.2f) ', ...
%         'in %0.0f trials. \n'], ...
%         ps_cv_start, cv_se(ps_cv_start, stim_trials(1)), ...
%         ps_skew_start, stim_skew_ses(1), ps_kurt_start, stim_kurt_ses(1), ...
%         stim_trials(1))
%     fprintf([...
%         'start-of-trial: CGS  induced increases in cv (%0.2f +/- %0.2f), '...
%         'skew (%0.2f +/- %0.2f) and kurtosis (%0.2f +/- %0.2f) ', ...
%         'in %0.0f trials. \n'], ...
%         gs_cv_start, cv_se(gs_cv_start, stim_trials(2)), ...
%         gs_skew_start, stim_skew_ses(2), gs_kurt_start, stim_kurt_ses(2), ...
%         stim_trials(2))
%     fprintf([...
%         'start-of-trial: CTRL induced increases in cv (%0.2f +/- %0.2f), '...
%         'skew (%0.2f +/- %0.2f) and kurtosis (%0.2f +/- %0.2f) ', ...
%         'in %0.0f trials. \n'], ...
%         ctrl_cv_start, cv_se(ctrl_cv_start, stim_trials(3)), ...
%         ctrl_skew_start, stim_skew_ses(3), ctrl_kurt_start, stim_kurt_ses(3), ...
%         stim_trials(3))
%     fprintf([...
%         'start-of-trial: AGS  induced increases in cv (%0.2f +/- %0.2f), '...
%         'skew (%0.2f +/- %0.2f) and kurtosis (%0.2f +/- %0.2f) ', ...
%         'in %0.0f trials. \n'], ...
%         an_cv_start, cv_se(an_cv_start, stim_trials(4)), ...
%         an_skew_start, stim_skew_ses(4), an_kurt_start, stim_kurt_ses(4), ...
%         stim_trials(4))
%     
%     % End of Trial
%     fprintf([...
%         'end-of-trial  : PS   induced increases in cv (%0.2f +/- %0.2f), '...
%         'skew (%0.2f +/- %0.2f) and kurtosis (%0.2f +/- %0.2f) ', ...
%         'in %0.0f trials. \n'], ...
%         ps_cv_end, cv_se(ps_cv_end, stim_trials(1)), ...
%         ps_skew_end, stim_skew_ses(1), ps_kurt_end, stim_kurt_ses(1), ...
%         stim_trials(1))
%     fprintf([...
%         'end-of-trial  : CGS  induced increases in cv (%0.2f +/- %0.2f), '...
%         'skew (%0.2f +/- %0.2f) and kurtosis (%0.2f +/- %0.2f) ', ...
%         'in %0.0f trials. \n'], ...
%         gs_cv_end, cv_se(gs_cv_end, stim_trials(2)), ...
%         gs_skew_end, stim_skew_ses(2), gs_kurt_end, stim_kurt_ses(2), ...
%         stim_trials(2))
%     fprintf([...
%         'end-of-trial  : CTRL induced increases in cv (%0.2f +/- %0.2f), '...
%         'skew (%0.2f +/- %0.2f) and kurtosis (%0.2f +/- %0.2f) ', ...
%         'in %0.0f trials. \n'], ...
%         ctrl_cv_end, cv_se(ctrl_cv_end, stim_trials(3)), ...
%         ctrl_skew_end, stim_skew_ses(3), ctrl_kurt_end, stim_kurt_ses(3), ...
%         stim_trials(3))
%     fprintf([...
%         'end-of-trial  : AGS  induced increases in cv (%0.2f +/- %0.2f), '...
%         'skew (%0.2f +/- %0.2f) and kurtosis (%0.2f +/- %0.2f) ', ...
%         'in %0.0f trials. \n'], ...
%         an_cv_end, cv_se(an_cv_end, stim_trials(4)), ...
%         an_skew_end, stim_skew_ses(4), an_kurt_end, stim_kurt_ses(4), ...
%         stim_trials(4))
    
    ball_rs = ball_rs(1:num_affected);
    ps_down = ps_frs(:, 1:num_affected);
    gs_down = gs_frs(:, 1:num_affected);
    ctrl_down = ctrl_frs(:, 1:num_affected);
    an_down = an_frs(:, 1:num_affected);
    
    % shift data so it is discontinuous
    if plot_name == "p1_wins"
        low_fr_lim = 50;
    else
        low_fr_lim = 15;
    end
    high_fr_lim = 200;
    frac_high = 1/3;
    ps_down(ps_down>=low_fr_lim) = low_fr_lim + (ps_down(ps_down>=low_fr_lim)-low_fr_lim) ./ ...
                                    (high_fr_lim - low_fr_lim) .* (frac_high*low_fr_lim);
    gs_down(gs_down>=low_fr_lim) = low_fr_lim + (gs_down(gs_down>=low_fr_lim)-low_fr_lim) ./ ...
                                    (high_fr_lim - low_fr_lim) .* (frac_high*low_fr_lim); 
    ctrl_down(ctrl_down>=low_fr_lim) = low_fr_lim + (ctrl_down(ctrl_down>=low_fr_lim)-low_fr_lim) ./ ...
                                    (high_fr_lim - low_fr_lim) .* (frac_high*low_fr_lim);
    an_down(an_down>=low_fr_lim) = low_fr_lim + (an_down(an_down>=low_fr_lim)-low_fr_lim) ./ ...
                                    (high_fr_lim - low_fr_lim) .* (frac_high*low_fr_lim);
    
    xlims = [0.5, 3.5]; %omit first 0.5s to avoid 0FR artifact 
    figure;
    hold on
    plot(t_down, ps_cv, 'Color', default_colors(7, :))
    plot(t_down, gs_cv, 'Color', default_colors(5, :))
    plot(t_down, an_cv, 'Color', default_colors(6, :))
    plot(t_down, ctrl_cv, 'k-')
    legend(["PS", "CGS", "AGS", "No Stim"])
    xlabel("Time (s)")
    xlim(xlims)
    ylim([0, 5])
    ylabel("Coefficient of Variation (unitless)")
    hold off
    
    figure;
    hold on
    plot(t_down, ps_skew, 'Color', default_colors(7, :))
    plot(t_down, gs_skew, 'Color', default_colors(5, :))
    plot(t_down, an_skew, 'Color', default_colors(6, :))
    plot(t_down, ctrl_skew, 'k')
    legend(["PS", "CGS", "AGS", "No Stim"])
    xlabel("Time (s)")
    xlim(xlims)
    ylim([-2, 10])
    ylabel("Skewness (unitless)")
    hold off
    
    figure;
    hold on
    plot(t_down, ps_kurt, 'Color', default_colors(7, : ))
    plot(t_down, gs_kurt, 'Color', default_colors(5, :))
    plot(t_down, an_kurt, 'Color', default_colors(6, :))
    plot(t_down, ctrl_kurt, 'k')
    legend(["PS", "CGS", "AGS", "No Stim"])
    xlabel("Time (s)")
    xlim(xlims)
    ylim([0, 100])
    ylabel("Kurtosis (unitless)")
    hold off
    
    near_mask = ball_rs<=500e-6;
    far_mask = logical(ones(size(ball_rs)));
    near_ticks = 0:100:500;
    far_ticks = 0:500:2000;
    near_map = pink;
    far_map = hot;
    near_lims = [0, 150];
    far_lims = [0, 50];
    clim = low_fr_lim + low_fr_lim * frac_high;
    %low_ticks = [0, 10, 20, 30, 40, 50];
    low_ticks = [0, 5, 10, 15];
    high_ticks = [50, 100, 150, 200];
    cticks = [low_ticks, low_fr_lim + (high_ticks-low_fr_lim)./(high_fr_lim-low_fr_lim).*low_fr_lim.*frac_high];
    cticklabels = [low_ticks, high_ticks];
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, ps_down(:, far_mask)');
    colormap(far_map)
    caxis([0, clim])
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    cbar.Ticks = cticks;
    cbar.TickLabels = cticklabels;
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("Pulse")
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, gs_down(:, far_mask)');
    colormap(far_map)
    caxis([0, clim])
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    cbar.Ticks = cticks;
    cbar.TickLabels = cticklabels;
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("Galvanic")
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, ctrl_down(:, far_mask)');
    colormap(far_map)
    caxis([0, clim])
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    cbar.Ticks = cticks;
    cbar.TickLabels = cticklabels;
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("Control")
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, an_down(:, far_mask)');
    colormap(far_map)
    caxis([0, clim])
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    cbar.Ticks = cticks;
    cbar.TickLabels = cticklabels;
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("Anodic")
    
    % Near vs Far version
    %{
    figure;
    s = pcolor(t_down, ball_rs(near_mask)*1e6, ps_down(:, near_mask)');
    colormap(near_map)
    caxis([0, 150])
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(near_ticks)
    title("Pulse Near")
    
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, ps_down(:, far_mask)');
    colormap(far_map)
    caxis(far_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("Pulse Far")
    
    figure;
    s = pcolor(t_down, ball_rs(near_mask)*1e6, gs_down(:, near_mask)');
    colormap(near_map)
    caxis(near_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(near_ticks)
    title("CGS Near")
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, ps_down(:, far_mask)');
    colormap(near_map)
    caxis(near_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("Pulse Near Map")
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, gs_down(:, far_mask)');
    colormap(far_map)
    caxis(far_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("CGS Far")
    
    figure;
    s = pcolor(t_down, ball_rs(near_mask)*1e6, ctrl_down(:, near_mask)');
    colormap(near_map)
    caxis(near_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(near_ticks)
    title("Control Near")
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, ctrl_down(:, far_mask)');
    colormap(far_map)
    caxis(far_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("Control Far")
    
    figure;
    s = pcolor(t_down, ball_rs(near_mask)*1e6, an_down(:, near_mask)');
    colormap(near_map)
    caxis(near_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(near_ticks)
    title("AGS Near")
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, an_down(:, far_mask)');
    colormap(far_map)
    caxis(far_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("AGS Far")
    %}
end

function neuron_frs = recspikes2neuron_frs(recspikes, win_size, avg_win_size, t, N)
    dt = t(2) - t(1);
    spikes = zeros(length(t), N);
    for nn = 1:N
        for spikeidx = recspikes(int2str(nn))
            spikes(spikeidx, nn) = 1;
        end
    end
    w = ones(floor(win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, spikes) ./ dt;
    w = ones(floor(avg_win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, neuron_frs);
end

function standard_error = cv_se(cv, n)
    % Computes the Standard Error of coefficient of variation (CV) based on
    % # of trials n AND cv
    % See https://influentialpoints.com/Training/standard_error_of_coefficient_of_variation.htm for derivation
    standard_error = cv ./ sqrt(2*n);
end

function standard_error = skew_se(n)
    % Computes the Standard Error of skewness based on # of trials n
    % See https://www.graphpad.com/support/faqid/1577/ for derivation
    standard_error = sqrt( 6*n.*(n-1) ./ ( (n-2).*(n+1).*(n+3) ) );
end

function standard_error = kurt_se(n)
    % Computes the Standard Error of kurtosis based on # of trials n
    % See https://www.youtube.com/watch?v=JEUYR_LjjqU for derivation
    standard_error = 2 * skew_se(n) .* sqrt(  (n.^2 - 1)./( (n-3).*(n+5) )  );
end