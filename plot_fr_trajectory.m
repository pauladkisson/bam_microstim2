%%% Paul Adkisson
%%% 12/6/2022
%%% Plot Mean Firing Rate Trajectories for each Stimulation Condition
function plot_fr_trajectory(sim_name, pulse_amps, stim_amps, t, t_cut, t_task, ...
    ex_c, pulse_coherences, galvanic_coherences, control_coherences, anodic_coherences, ...
    default_colors, start_trial, end_trial, num_trials, N, p, f, N_E, ...
    start_thresh, stop_thresh, plot_name)
    if plot_name == "p1_wins"
        win_num = 1;
    else
        win_num = 2;
    end
    dt = t(2) - t(1);
    num_group = floor(f*N_E);
    stim_frs = zeros(length(stim_amps), num_trials, length(t));
    aligned_t = dt-1:dt:1-dt;
    aligned_frs = zeros(length(stim_amps), num_trials, length(aligned_t));
    stim_slopes = zeros(length(stim_amps), num_trials);
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
                aligned_frs(j, relative_trial, :) = NaN;
                stim_slopes(j, relative_trial) = NaN;
                continue %skip trials where P1 doesn't win/lose or decision takes too long
            end
            load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                   "recspikes")
            [pop_frs, ~] = recspikes2popfrs(recspikes, t, N, dt, p, f, N_E);
            stim_frs(j, relative_trial, :) = pop_frs(:, 1);
            aligned_mask = t>=dec_time+t_task-1+dt/2 & ...
                           t<dec_time+t_task+1-dt/2;
            aligned_frs(j, relative_trial, :) = pop_frs(aligned_mask, 1);
            if plot_name == "p1_wins"
                preidx = find(pop_frs(:, 1)>=start_thresh, 1);
                postidx = find(pop_frs(:, 1)>=stop_thresh, 1);
                slope_t = t(preidx:postidx);
                slope = (stop_thresh-start_thresh) / (t(postidx)-t(preidx));
                coeffs = [slope, pop_frs(preidx, 1) - slope*t(preidx)];
                %coeffs = polyfit(slope_t, pop_frs(preidx:postidx, 1), 1);
                stim_slopes(j, relative_trial) = coeffs(1);
                
                %plot to debug outliers
                debug = coeffs(1) > 100 | coeffs(1) < 30;
                if debug
                    figure;
                    slope_y = coeffs(2) + coeffs(1)*slope_t;
                    slope_y = slope_y - (slope_y(2) - slope_y(1))/2; %center
                    hold on;
                    plot(t, pop_frs(:, 1))
                    plot(slope_t, slope_y, "r--")
                    title(sprintf("j=%0.0f", j))
                end                     
            end
        end
    end
    pulse_frs = reshape(stim_frs(1, :, :), [num_trials, length(t)]);
    galvanic_frs = reshape(stim_frs(2, :, :), [num_trials, length(t)]);
    control_frs = reshape(stim_frs(3, :, :), [num_trials, length(t)]);
    anodic_frs = reshape(stim_frs(4, :, :), [num_trials, length(t)]);
    pulse_trialmean = mean(pulse_frs, 1, 'omitnan');
    galvanic_trialmean = mean(galvanic_frs, 1, 'omitnan');
    control_trialmean = mean(control_frs, 1, 'omitnan');
    anodic_trialmean = mean(anodic_frs, 1, 'omitnan');
    
    figure;
    set(gca, 'Fontsize', 18)
    hold on
    plot(t, pulse_trialmean, 'Color', default_colors(7, :), 'Linewidth', 2)
    plot(t, galvanic_trialmean, 'Color', default_colors(5, :), 'Linewidth', 2)
    plot(t, control_trialmean, "k", 'Linewidth', 2)
    plot(t, anodic_trialmean, 'Color', default_colors(6, :), 'Linewidth', 2)
    hold off
    xlabel("Time (s)")
    ylabel("P1 Firing Rate (spk/s)")
    if plot_name == "p1_wins"
        title("P1 Wins")
    else
        title("P1 Loses")
    end
    
    % Slope around decision threshold
    pulse_slopes = reshape(stim_slopes(1, :, :), [num_trials, 1]);
    galvanic_slopes = reshape(stim_slopes(2, :, :), [num_trials, 1]);
    control_slopes = reshape(stim_slopes(3, :, :), [num_trials, 1]);
    anodic_slopes = reshape(stim_slopes(4, :, :), [num_trials, 1]);
    mean_ctrl = mean(control_slopes, 'omitnan');
    norm_ps = pulse_slopes - mean_ctrl;
    norm_gs = galvanic_slopes - mean_ctrl;
    norm_ctrl = control_slopes - mean_ctrl;
    norm_an = anodic_slopes - mean_ctrl;
    ps_quantiles = quantile(norm_ps, [0.25, 0.5, 0.75]);
    gs_quantiles = quantile(norm_gs, [0.25, 0.5, 0.75]);
    an_quantiles = quantile(norm_an, [0.25, 0.5, 0.75]);
    ctrl_quantiles = quantile(norm_ctrl, [0.25, 0.5, 0.75]);
    
    figure;
    set(gca, 'fontsize', 18)
    hold on
    colors = [default_colors(5, :); default_colors(6, :); default_colors(7, :); [0, 0, 0]];
    boxplot([norm_gs, norm_an, norm_ps, norm_ctrl], 'PlotStyle', 'traditional', ...
        'Colors', colors, 'Symbol', ".")
    hold off
    xticks([1, 2, 3, 4])
    xticklabels(["Galvanic", "Anodic", "Pulsatile", "Control"])
    ylim([-40, 120])
    ylabel("Change in Firing Rate Slope (spk/s^2)")
    title("Recurrent Excitation Metric: P1 Wins")
    
    % Statistics
    disp("SLOPES")
    [p_median, ~, stats] = kruskalwallis([norm_ps, norm_gs, norm_ctrl, norm_an]);
    fprintf([...
        "Stimulation induces significantly different median slopes (p=%0.1e). \n"], ...
        p_median)
    c = multcompare(stats);
    p_ps_cgs_median = c(1, end);
    p_cgs_ctrl_median = c(4, end);
    p_ps_ctrl_median = c(2, end);
    p_cgs_ags_median = c(5, end);
    fprintf([...
        'CGS (%0.2f, %0.2f, %0.2f) has a higher median slope than ', ...
        'control (%0.2f, %0.2f, %0.2f) (p=%0.2f). \n'], ...
        gs_quantiles(1), gs_quantiles(2), gs_quantiles(3), ...
        ctrl_quantiles(1), ctrl_quantiles(2), ctrl_quantiles(3), ...
        p_cgs_ctrl_median)
    fprintf([...
        'CGS (%0.2f, %0.2f, %0.2f) has a higher median slope than ', ...
        'pulse (%0.2f, %0.2f, %0.2f) (p=%0.2f). \n'], ...
        gs_quantiles(1), gs_quantiles(2), gs_quantiles(3), ...
        ps_quantiles(1), ps_quantiles(2), ps_quantiles(3), ...
        p_ps_cgs_median)
    fprintf([...
        'PS (%0.2f, %0.2f, %0.2f) has a lower median slope than ', ...
        'control (%0.2f, %0.2f, %0.2f) (p=%0.2f). \n'], ...
        ps_quantiles(1), ps_quantiles(2), ps_quantiles(3), ...
        ctrl_quantiles(1), ctrl_quantiles(2), ctrl_quantiles(3), ...
        p_ps_ctrl_median)
    fprintf([...
        'CGS (%0.2f, %0.2f, %0.2f) has a higher median slope than ', ...
        'AGS (%0.2f, %0.2f, %0.2f) (p=%0.2f). \n'], ...
        gs_quantiles(1), gs_quantiles(2), gs_quantiles(3), ...
        an_quantiles(1), an_quantiles(2), an_quantiles(3), ...
        p_cgs_ags_median)
       
    % Max FRS
    pulse_max_frs = max(pulse_frs, [], 2);
    galvanic_max_frs = max(galvanic_frs, [], 2);
    control_max_frs = max(control_frs, [], 2);
    anodic_max_frs = max(anodic_frs, [], 2);
    mean_ctrl = mean(control_max_frs, 'omitnan');
    norm_ps = pulse_max_frs - mean_ctrl;
    norm_gs = galvanic_max_frs - mean_ctrl;
    norm_an = anodic_max_frs - mean_ctrl;
    norm_ctrl = control_max_frs - mean_ctrl;
    ps_quantiles = quantile(norm_ps, [0.25, 0.5, 0.75]);
    gs_quantiles = quantile(norm_gs, [0.25, 0.5, 0.75]);
    an_quantiles = quantile(norm_an, [0.25, 0.5, 0.75]);
    ctrl_quantiles = quantile(norm_ctrl, [0.25, 0.5, 0.75]);
        
    figure;
    set(gca, 'fontsize', 18)
    hold on
    boxplot([norm_gs, norm_an, norm_ps, norm_ctrl], 'PlotStyle', 'traditional', ...
        'Colors', colors, 'Symbol', ".")
    hold off
    xticks([1, 2, 3, 4])
    xticklabels(["Galvanic", "Anodic", "Pulsatile", "Control"])
    ylim([-5, 15])
    ylabel("Change in Maximum Firing Rate (spk/s)")
    title("Recurrent Excitation Metric: P1 Loses")
    
    % Statistics
    [p_median, ~, stats] = kruskalwallis([norm_ps, norm_gs, norm_ctrl, norm_an]);
    fprintf([...
        'Stimulation induces significantly different median Max FR (p=%0.1e). \n'], ...
        p_median)
    c = multcompare(stats);
    p_ps_cgs_median = c(1, end);
    p_cgs_ctrl_median = c(4, end);
    p_ags_ctrl_median = c(6, end);
    p_ps_ctrl_median = c(2, end);
    fprintf([...
        'CGS (%0.2f, %0.2f, %0.2f) has a higher median maximum FR than ', ...
        'control (%0.2f, %0.2f, %0.2f) (p=%0.1e). \n'], ...
        gs_quantiles(1), gs_quantiles(2), gs_quantiles(3), ...
        ctrl_quantiles(1), ctrl_quantiles(2), ctrl_quantiles(3), ...
        p_cgs_ctrl_median)
    fprintf([...
        'PS (%0.2f, %0.2f, %0.2f) has a higher median maximum FR than ', ...
        'control (%0.2f, %0.2f, %0.2f) (p=%0.1e). \n'], ...
        ps_quantiles(1), ps_quantiles(2), ps_quantiles(3), ...
        ctrl_quantiles(1), ctrl_quantiles(2), ctrl_quantiles(3), ...
        p_ps_ctrl_median)
    fprintf([...
        'AGS (%0.2f, %0.2f, %0.2f) has a lower median maximum FR than ', ...
        'control (%0.2f, %0.2f, %0.2f) (p=%0.1e). \n'], ...
        an_quantiles(1), an_quantiles(2), an_quantiles(3), ...
        ctrl_quantiles(1), ctrl_quantiles(2), ctrl_quantiles(3), ...
        p_ags_ctrl_median)
    fprintf('PS and CGS have equivalent median maximum FRs (p=%0.2f). \n', p_ps_cgs_median)
end