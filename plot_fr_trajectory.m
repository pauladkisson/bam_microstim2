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
                %slope = (stop_thresh-start_thresh) / (t(postidx)-t(preidx));
                %coeffs = [slope, pop_frs(preidx, 1) - slope*t(preidx)];
                coeffs = polyfit(slope_t, pop_frs(preidx:postidx, 1), 1);
                stim_slopes(j, relative_trial) = coeffs(1);
                
                %plot to debug outliers
                debug = coeffs(1) > 100; 
                if debug
                    figure;
                    slope_y = coeffs(2) + coeffs(1)*slope_t;
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
    
    pulse_slopes = reshape(stim_slopes(1, :, :), [num_trials, 1]);
    galvanic_slopes = reshape(stim_slopes(2, :, :), [num_trials, 1]);
    control_slopes = reshape(stim_slopes(3, :, :), [num_trials, 1]);
    anodic_slopes = reshape(stim_slopes(4, :, :), [num_trials, 1]);
    mean_ctrl = mean(control_slopes, 'omitnan');
    norm_ps = pulse_slopes - mean_ctrl;
    norm_gs = galvanic_slopes - mean_ctrl;
    norm_an = anodic_slopes - mean_ctrl;
    mean_slopes = [mean(norm_gs, 'omitnan'), mean(norm_an, 'omitnan'), mean(norm_ps, 'omitnan')];
    std_slopes = [std(norm_gs, [], 'omitnan'), std(norm_an, [], 'omitnan'), std(norm_ps, [], 'omitnan')];
    stim_trials = [sum(~isnan(norm_gs)), sum(~isnan(norm_an)), sum(~isnan(norm_ps))];
    
    % Statistics
    disp("SLOPES")
    stim_slopes = [pulse_slopes, galvanic_slopes, control_slopes, anodic_slopes];
    [~, ~, stats] = anova1(stim_slopes, [], 'off');
    results = multcompare(stats, 'Display', 'off');
    p_vals = results(:, end);
    p_ps = p_vals(2);
    p_gs = p_vals(4);
    p_ags = p_vals(end);
    fprintf([...
        'CGS (p=%0.2f) increased the slope by %0.2f +/- %0.2f ', ...
        'relative to control. \n'], ...
        p_gs, mean_slopes(1), std_slopes(1)/sqrt(stim_trials(1)))
    fprintf([...
        'PS (p=%0.2f) and AGS (p=%0.2f) decreased the slope ', ...
        'by %0.2f +/- %0.2f and %0.2f +/- %0.2f respectively ', ...
        'relative to control. \n'], p_ps, p_ags, mean_slopes(3), ...
         std_slopes(3)/sqrt(stim_trials(3)), mean_slopes(2), ...
         std_slopes(2)/sqrt(stim_trials(2)))
    
    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(mean_slopes);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(6, :); default_colors(7, :)];
    %b.CData = [default_colors(7, :); default_colors(5, :); [0, 0, 0]; default_colors(6, :)];
    x = [1, 2, 3];
    errorbar(x, mean_slopes, std_slopes, 'k.', 'Linewidth', 20, 'Capsize', 0)
    hold off
    xticks([1, 2, 3])
    xticklabels(["Galvanic", "Anodic", "Pulsatile"])
    ylabel("Change in Firing Rate Slope (spk/s^2)")
    title("Recurrent Excitation Metric")
    
    pulse_max_frs = max(pulse_frs, [], 2);
    galvanic_max_frs = max(galvanic_frs, [], 2);
    control_max_frs = max(control_frs, [], 2);
    anodic_max_frs = max(anodic_frs, [], 2);
    mean_ctrl = mean(control_max_frs, 'omitnan');
    norm_ps = pulse_max_frs - mean_ctrl;
    norm_gs = galvanic_max_frs - mean_ctrl;
    norm_an = anodic_max_frs - mean_ctrl;
    mean_max_frs = [mean(norm_gs, 'omitnan'), mean(norm_an, 'omitnan'), mean(norm_ps, 'omitnan')];
    std_max_frs = [std(norm_gs, [], 'omitnan'), std(norm_an, [], 'omitnan'), std(norm_ps, [], 'omitnan')];
    
    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(mean_max_frs);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(6, :); default_colors(7, :)];
    x = [1, 2, 3];
    errorbar(x, mean_max_frs, std_max_frs, 'k.', 'Linewidth', 20, 'Capsize', 0)
    hold off
    xticks([1, 2, 3])
    xticklabels(["Galvanic", "Anodic", "Pulsatile"])
    ylabel("Change in Maximum Firing Rate (spk/s)")
    title("Recurrent Excitation Metric: P1 Loses")
    
    % Statistics
    disp("MAX FRS")
    [~, p_ps_cgs] = ttest2(norm_ps, norm_gs);
    fprintf([...
        'Maximum P1 FR was increased by PS (%0.2f +/- %0.2f) and CGS ', ...
        '(%0.2f +/- %0.2f) but decreased by AGS (%0.2f +/- %0.2f) ', ...
        'relative to control. \n'], ...
        mean_max_frs(3), std_max_frs(3)/sqrt(stim_trials(3)), ...
        mean_max_frs(1), std_max_frs(1)/sqrt(stim_trials(1)), ...
        mean_max_frs(2), std_max_frs(2)/sqrt(stim_trials(2)))
    fprintf('PS and CGS induced statistically equivalent increases (p=%0.2f)', ...
         p_ps_cgs)
    
    ps_aligned = reshape(aligned_frs(1, :, :), [num_trials, length(aligned_t)]);
    ps_aligned_mean = mean(ps_aligned, 1, 'omitnan');
    gs_aligned = reshape(aligned_frs(2, :, :), [num_trials, length(aligned_t)]);
    gs_aligned_mean = mean(gs_aligned, 1, 'omitnan');
    ctrl_aligned = reshape(aligned_frs(3, :, :), [num_trials, length(aligned_t)]);
    ctrl_aligned_mean = mean(ctrl_aligned, 1, 'omitnan');
    an_aligned = reshape(aligned_frs(4, :, :), [num_trials, length(aligned_t)]);
    an_aligned_mean = mean(an_aligned, 1, 'omitnan');
    
    figure;
    hold on
    plot(aligned_t, ps_aligned, 'Color', default_colors(7, :))
    plot(aligned_t, ps_aligned_mean, 'Color', default_colors(7, :), 'Linewidth', 2)
    plot(aligned_t, gs_aligned, 'Color', default_colors(5, :))
    plot(aligned_t, gs_aligned_mean, 'Color', default_colors(5, :), 'Linewidth', 2)
    plot(aligned_t, an_aligned, 'Color', default_colors(6, :))
    plot(aligned_t, an_aligned_mean, 'Color', default_colors(6, :), 'Linewidth', 2)
    plot(aligned_t, ctrl_aligned, 'k-')
    plot(aligned_t, ctrl_aligned_mean, 'k-', 'Linewidth', 2)
    hold off
end