%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Decision Metrics Accuracy and Decision Time
function plot_decisions(sim_name, pulse_amps, stim_amps, t_cut, ...
                default_colors, num_batch, num_trials, pulse_coherences, ...
                galvanic_coherences, control_coherences, anodic_coherences)
    %Accuracy
    num_amps = length(stim_amps);
    stim_coeffs = zeros(num_amps, num_batch, 2);
    c = -1:0.001:1;
    stim_ws = zeros(num_amps, num_batch, length(c));
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j <= length(pulse_amps);
        if pulse
            datapath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                [sim_name, stim_amp*1e6]);
            stim_coherences = pulse_coherences;
        else
            datapath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e6]);
            if stim_amp < 0 %cathodic gs
                stim_coherences = galvanic_coherences;
            elseif stim_amp == 0 %Control
                stim_coherences = control_coherences;
            else %anodic gs
                stim_coherences = anodic_coherences;
            end
        end
        load(strcat(datapath, "/decisions.mat"), "avg_acc", ...
            "coeffs", "batch_coeffs", "batch_acc", "percent_nodec");
        if pulse
            pulse_acc = avg_acc;
            pulse_acc_sem = std(batch_acc, [], 2) ./ sqrt(num_batch);
            pulse_coeffs = coeffs;
            pulse_nodec = percent_nodec;
        elseif stim_amp < 0
            galvanic_acc = avg_acc;
            galvanic_acc_sem = std(batch_acc, [], 2) ./ sqrt(num_batch);
            galvanic_coeffs = coeffs;
            galvanic_nodec = percent_nodec;
        elseif stim_amp == 0
            ctrl_acc = avg_acc;
            ctrl_acc_sem = std(batch_acc, [], 2) ./ sqrt(num_batch);
            ctrl_coeffs = coeffs;
            ctrl_nodec = percent_nodec;
        else %anodic
            anodic_acc = avg_acc;
            anodic_acc_sem = std(batch_acc, [], 2) ./ sqrt(num_batch);
            anodic_coeffs = coeffs;
            anodic_nodec = percent_nodec;
        end
        stim_coeffs(j, :, :) = batch_coeffs;
    end
    pulse_w = logistic_acc(pulse_coeffs, c);
    galvanic_w = logistic_acc(galvanic_coeffs , c);
    control_w = logistic_acc(ctrl_coeffs , c);
    anodic_w = logistic_acc(anodic_coeffs , c);
    
    figure;
    set(gca, 'fontsize', 18);
    hold on
    errorbar(control_coherences, ctrl_acc, ctrl_acc_sem, 'k.', 'MarkerSize', 20)
    errorbar(pulse_coherences, pulse_acc, pulse_acc_sem, '.', ...
        'Color', default_colors(7, :), 'MarkerSize', 20)
    errorbar(galvanic_coherences, galvanic_acc, galvanic_acc_sem, '.', ...
        'Color', default_colors(5, :), 'MarkerSize', 20)
    errorbar(anodic_coherences, anodic_acc, anodic_acc_sem, '.', ...
         'Color', default_colors(6, :), 'MarkerSize', 20)
    plot(c, control_w, "k")
    plot(c, galvanic_w, 'Color', default_colors(5, :))
    plot(c, anodic_w, 'Color', default_colors(6, :))
    plot(c, pulse_w, 'Color', default_colors(7, :))
    hold off
    xlabel("Coherence (%)")
    ylabel("% of trials P1 wins")
    legend("Pulsatile", "Galvanic", "Control", "Anodic")
    ylim([0, 1])
    
    % Compute Psychometric Parameters
    a = reshape(stim_coeffs(:, :, 1), [length(stim_amps), num_batch]);
    b = reshape(stim_coeffs(:, :, 2), [length(stim_amps), num_batch]);
    bias = -a ./ b;
    sensitivity = b ./ 4;

    % Bias Stats
    avg_ps_bias = mean(bias(1, :));
    sem_ps_bias = std(bias(1, :)) ./ sqrt(num_batch);
    avg_gs_bias = mean(bias(2, :));
    sem_gs_bias = std(bias(2, :)) ./ sqrt(num_batch);
    avg_an_bias = mean(bias(4, :));
    sem_an_bias = std(bias(4, :)) ./ sqrt(num_batch);
    [~, p_cgs_ags] = ttest2(abs(bias(2, :)), abs(bias(4, :)));
    [~, p_ps_cgs] = ttest2(bias(1, :), bias(2, :));
    [~, p_ctrl] = ttest(bias(3, :));
    disp("BIAS")
    fprintf([...
        'PS and CGS shifted psychometric curve by %0.1f +/- %0.1f%%', ...
        ' and %0.1f +/- %0.1f%% respectively. \n'], ...
        avg_ps_bias*100, sem_ps_bias*100, avg_gs_bias*100, sem_gs_bias*100)
    fprintf('AGS shifted the psychometric curve by %0.1f +/- %0.1f%%. \n', ...
            avg_an_bias*100, sem_an_bias*100)
    fprintf([...
        'Control did not significantly bias the psychometric curve ', ...
        '(p = %0.2f). \n'], p_ctrl)
    fprintf([...
        'AGS shifted the psychometric curve significantly less than CGS ', ...
        '(p=%0.2e). \n'], p_cgs_ags)
    fprintf([...
        'PS and CGS produced statistically identical bias (p=%0.2f). \n'], ...
       p_ps_cgs)

   % Sensitivity Stats
    ps_sens = reshape(sensitivity(1, :), [num_batch, 1]);
    gs_sens = reshape(sensitivity(2, :), [num_batch, 1]);
    ctrl_sens = reshape(sensitivity(3, :), [num_batch, 1]);
    an_sens = reshape(sensitivity(4, :), [num_batch, 1]);
    mean_ctrl_sens = mean(ctrl_sens);
    norm_ps_sens = ps_sens - mean_ctrl_sens;
    norm_gs_sens = gs_sens - mean_ctrl_sens;
    norm_an_sens = an_sens - mean_ctrl_sens;
    
    avg_norm_ps_sens = mean(norm_ps_sens);
    avg_norm_gs_sens = mean(norm_gs_sens);
    avg_norm_an_sens = mean(norm_an_sens);
    sem_norm_ps_sens = std(norm_ps_sens) ./ sqrt(num_batch);
    sem_norm_gs_sens = std(norm_gs_sens) ./ sqrt(num_batch);
    sem_norm_an_sens = std(norm_an_sens) ./ sqrt(num_batch);
    [~, ~, stats] = anova1(sensitivity', [], 'off');
    results = multcompare(stats, 'Display', 'off');
    p_vals = results(:, end);
    p_ps = p_vals(2);
    p_gs = p_vals(4);
    p_ags = p_vals(end);
    
    disp("SENSITIVITY")
    fprintf([...
        'PS (p=%0.2f) and CGS (p=%0.2f) decreased the sensitivity ', ...
        'by %0.2f +/- %0.2f and %0.2f +/- %0.2f respectively ', ...
        'relative to control. \n'], p_ps, p_gs, avg_norm_ps_sens, ...
        sem_norm_ps_sens, avg_norm_gs_sens, sem_norm_gs_sens)
    
    fprintf([...
        'AGS (p=%0.2f) increased the sensitivity by %0.2f +/- %0.2f ', ...
        'relative to control. \n'], ...
        p_ags, avg_norm_an_sens, sem_norm_an_sens)

    
    %Decision Times
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j <= length(pulse_amps);
        if pulse
            datapath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                [sim_name, stim_amp*1e6]);
        else
            datapath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e6]);
        end
        load(strcat(datapath, "/decisions.mat"), "decision_times");
        if pulse
            pulse_dt = decision_times;
            pulse_latedec = sum(decision_times > t_cut, 1) / num_trials;
        elseif stim_amp < 0 %cathodic GS
            galvanic_dt = decision_times;
            galvanic_latedec = sum(decision_times > t_cut, 1) / num_trials;
        elseif stim_amp == 0 %control
            control_dt = decision_times;
            ctrl_latedec = sum(decision_times > t_cut, 1) / num_trials;
        else %anoidc GS
            anodic_dt = decision_times;
            anodic_latedec = sum(decision_times > t_cut, 1) / num_trials;
        end
    end
    trialmean_pulse = mean(pulse_dt, 1, 'omitnan');
    trialstd_pulse = std(pulse_dt, [], 1, 'omitnan');
    trialmean_galvanic = mean(galvanic_dt, 1, 'omitnan');
    trialstd_galvanic = std(galvanic_dt, [], 1, 'omitnan');
    trialmean_ctrl = mean(control_dt, 1, 'omitnan');
    trialstd_ctrl = std(control_dt, [], 1, 'omitnan');
    trialmean_anodic = mean(anodic_dt, 1, 'omitnan');
    trialstd_anodic = std(anodic_dt, 1, 'omitnan');
    figure;
    set(gca, 'fontsize', 18);
    hold on
    errorbar(pulse_coherences, trialmean_pulse, trialstd_pulse/sqrt(num_trials), 'Color', default_colors(7, :))
    errorbar(galvanic_coherences, trialmean_galvanic, trialstd_galvanic/sqrt(num_trials), 'Color', default_colors(5, :))
    errorbar(control_coherences, trialmean_ctrl, trialstd_ctrl/sqrt(num_trials), 'k')
    errorbar(anodic_coherences, trialmean_anodic, trialstd_anodic/sqrt(num_trials), 'Color', default_colors(6, :))
    xticks([-1, -0.5, 0, 0.5, 1])
    ylim([0, 1.6])
%     scatter(control_coherences, ctrl_dt', 'k', 'filled')
%     scatter(pulse_coherences, pulse_dt', [], default_colors(7, :).*ones(length(pulse_dt), 3), 'filled')
%     scatter(galvanic_coherences, galvanic_dt', [], default_colors(5, :).*ones(length(galvanic_dt), 3), 'filled')
%     scatter(anodic_coherences, anodic_dt', [], default_colors(6, :).*ones(length(anodic_dt), 3), 'filled')
    hold off
    xlabel("Coherence (%)")
    ylabel("Decision Time (s)")
    legend("Pulsatile", "Galvanic", "Control", "Anodic")
    
    % DT stats
    [ps_peak_dt, ps_peak_coh] = max(pulse_dt, [], 2, 'omitnan');
    ps_peak_coh = pulse_coherences(ps_peak_coh);
    [gs_peak_dt, gs_peak_coh] = max(galvanic_dt, [], 2, 'omitnan');
    gs_peak_coh = galvanic_coherences(gs_peak_coh);
    [ctrl_peak_dt, ctrl_peak_coh] = max(control_dt, [], 2, 'omitnan');
    ctrl_peak_coh = control_coherences(ctrl_peak_coh); 
    [an_peak_dt, an_peak_coh] = max(anodic_dt, [], 2, 'omitnan');
    an_peak_coh = anodic_coherences(an_peak_coh);
    avg_ctrl_peak_dt = mean(ctrl_peak_dt);
    avg_ctrl_peak_coh = mean(ctrl_peak_coh);
    norm_ps_peak_dt = ps_peak_dt - avg_ctrl_peak_dt;
    norm_ps_peak_coh = ps_peak_coh - avg_ctrl_peak_coh;
    norm_gs_peak_dt = gs_peak_dt - avg_ctrl_peak_dt;
    norm_gs_peak_coh = gs_peak_coh - avg_ctrl_peak_coh;
    norm_an_peak_dt = an_peak_dt - avg_ctrl_peak_dt;
    norm_an_peak_coh = an_peak_coh - avg_ctrl_peak_coh;
    
    avg_norm_ps_peak_dt = mean(norm_ps_peak_dt);
    sem_norm_ps_peak_dt = std(norm_ps_peak_dt) ./ sqrt(length(ps_peak_dt));
    avg_norm_gs_peak_dt = mean(norm_gs_peak_dt);
    sem_norm_gs_peak_dt = std(norm_gs_peak_dt) ./ sqrt(length(gs_peak_dt));
    avg_norm_an_peak_dt = mean(norm_an_peak_dt);
    sem_norm_an_peak_dt = std(norm_an_peak_dt) ./ sqrt(length(an_peak_dt));
    avg_norm_ps_peak_coh = mean(norm_ps_peak_coh);
    sem_norm_ps_peak_coh = std(norm_ps_peak_coh) ./ sqrt(length(ps_peak_coh));
    avg_norm_gs_peak_coh = mean(norm_gs_peak_coh);
    sem_norm_gs_peak_coh = std(norm_gs_peak_coh) ./ sqrt(length(gs_peak_coh));
    avg_norm_an_peak_coh = mean(norm_an_peak_coh);
    sem_norm_an_peak_coh = std(norm_an_peak_coh) ./ sqrt(length(an_peak_coh));
    [~, p_cgs_ags_dt] = ttest2(abs(norm_gs_peak_dt), abs(norm_an_peak_dt));
    [~, p_cgs_ags_coh] = ttest2(abs(norm_gs_peak_coh), abs(norm_an_peak_coh));
    [~, p_ps_cgs_dt] = ttest2(norm_ps_peak_dt, norm_gs_peak_dt);
    [~, p_ps_cgs_coh] = ttest2(norm_ps_peak_coh, norm_gs_peak_coh);
    
    
    disp("DECISION TIME")
    fprintf([...
        'PS and CGS decreased peak decision time by %0.2f +/- %0.2fs and ', ...
        '%0.2f +/- %0.2fs respectively. \n'], avg_norm_ps_peak_dt, ...
        sem_norm_ps_peak_dt, avg_norm_gs_peak_dt, sem_norm_gs_peak_dt)
    fprintf([...
        'PS and CGS decreased peak coherence by %0.1f +/- %0.1f%% and ', ...
        '%0.1f +/- %0.1f%% respectively. \n'], avg_norm_ps_peak_coh*100, ...
        sem_norm_ps_peak_coh*100, avg_norm_gs_peak_coh*100, sem_norm_gs_peak_coh*100)
    fprintf([...
        'AGS increased peak decision time by %0.2f +/- %0.2f. \n'], ...
        avg_norm_an_peak_dt, sem_norm_an_peak_dt)
    fprintf([...
        'AGS increased peak coherence by %0.1f +/- %0.1f%%. \n'], ...
        avg_norm_an_peak_coh*100, sem_norm_an_peak_coh*100)
    fprintf([...
        'AGS shifted peak decision time (p=%0.2f) and coherence (%0.1e) ', ...
        'less than CGS. \n'], p_cgs_ags_dt, p_cgs_ags_coh)
    fprintf([...
        'PS and CGS induced statistically equivalent effects on ', ...
        'peak decision time (p=%0.2f) and peak coherence (p=%0.2f). \n'], ...
        p_ps_cgs_dt, p_ps_cgs_coh)
    
    %No-decisions
    figure;
    set(gca, 'fontsize', 18);
    hold on
    scatter(pulse_coherences, pulse_nodec*100, [], default_colors(7, :).*ones(length(pulse_nodec), 3), 'filled')
    scatter(galvanic_coherences, galvanic_nodec*100, [], default_colors(5, :).*ones(length(galvanic_nodec), 3), 'filled')
    scatter(anodic_coherences, anodic_nodec*100, [], default_colors(6, :).*ones(length(anodic_nodec), 3), 'filled')
    scatter(control_coherences, ctrl_nodec*100, 'k', 'filled') 
    hold off
    ylim([0, 100])
    xlabel("Coherence (%)")
    ylabel("% of trials No-Decision")
    
    %Late-decisions
    figure;
    set(gca, 'fontsize', 18);
    hold on
    scatter(pulse_coherences, pulse_latedec*100, [], default_colors(7, :).*ones(length(pulse_latedec), 3), 'filled')
    scatter(galvanic_coherences, galvanic_latedec*100, [], default_colors(5, :).*ones(length(galvanic_latedec), 3), 'filled')
    scatter(anodic_coherences, anodic_latedec*100, [], default_colors(6, :).*ones(length(anodic_latedec), 3), 'filled')
    scatter(control_coherences, ctrl_latedec*100, 'k', 'filled') 
    hold off
    ylim([0, 100])
    xlabel("Coherence (%)")
    ylabel("% of trials Late-Decision")
    title(sprintf("t_cut = %0.2f", t_cut))
end