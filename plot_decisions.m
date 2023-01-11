%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Decision Metrics Accuracy and Decision Time
function plot_decisions(sim_name, pulse_amps, stim_amps, default_colors, ...
                num_batch, num_trials, pulse_coherences, ...
                galvanic_coherences, control_coherences, anodic_coherences)
    %Accuracy
    num_amps = length(stim_amps);
    stim_coeffs = zeros(num_amps, num_batch, 2);
    c = -1:0.01:1;
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
        load(strcat(datapath, "/decisions.mat"), "avg_acc", "decisions", ...
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
        [nodec_trial, nodec_c] = find(~decisions);
        for nodec = 1:length(nodec_trial)
            fprintf("c=%0.3f, trial = %0.0f \n", ...
                [stim_coherences(nodec_c(nodec)), nodec_trial(nodec)])
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
    %{
    errorbar(pulse_coherences, pulse_acc, pulse_acc_sem, '.', ...
        'Color', default_colors(7, :), 'MarkerSize', 20)
    errorbar(galvanic_coherences, galvanic_acc, galvanic_acc_sem, '.', ...
        'Color', default_colors(5, :), 'MarkerSize', 20)
    errorbar(anodic_coherences, anodic_acc, anodic_acc_sem, '.', ...
        'Color', default_colors(6, :), 'MarkerSize', 20)
    %}
    plot(c, control_w, "k")
    %{ 
    plot(c, galvanic_w, 'Color', default_colors(5, :))
    plot(c, anodic_w, 'Color', default_colors(6, :))
    plot(c, pulse_w, 'Color', default_colors(7, :))
    hold off
    xlabel("Coherence (%)")
    ylabel("% of trials P1 wins")
    legend("Pulsatile", "Galvanic", "Control", "Anodic")
    ylim([0, 1])
    %}
    
    %Statistics
    a = reshape(stim_coeffs(:, 2:end, 1), [length(stim_amps), num_batch-1]);
    b = reshape(stim_coeffs(:, 2:end, 2), [length(stim_amps), num_batch-1]);
    bias = -a ./ b;
    sensitivity = b .* exp(-a) / (1 + exp(-a)).^2;

    avg_ps_bias = mean(bias(1, :))
    sem_ps_bias = std(bias(1, :)) ./ sqrt(num_batch)
    avg_gs_bias = mean(bias(2, :))
    sem_gs_bias = std(bias(2, :)) ./ sqrt(num_batch)
    avg_an_bias = mean(bias(4, :))
    sem_an_bias = std(bias(4, :)) ./ sqrt(num_batch)

    avg_ps_sens = mean(sensitivity(1, :))
    sem_ps_sens = std(sensitivity(1, :))./ sqrt(num_batch)
    avg_an_sens = mean(sensitivity(4, :))
    sem_an_sens = std(sensitivity(4, :)) ./ sqrt(num_batch)
    dc = c(2) - c(1);
    center_c0 = abs(c - avg_an_bias) < dc/2;
    center_c1 = abs(c - avg_an_bias + dc) < dc/2;
    num_an_sens = (anodic_w(center_c1) - anodic_w(center_c0)) / dc
    num_ctrl_sens = (control_w(c==dc) - control_w(c==0)) / dc
    
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
        elseif stim_amp < 0 %cathodic GS
            galvanic_dt = decision_times;
        elseif stim_amp == 0 %control
            ctrl_dt = decision_times;
        else %anoidc GS
            anodic_dt = decision_times;
        end
    end
    trialmean_pulse = mean(pulse_dt, 1, 'omitnan');
    trialstd_pulse = std(pulse_dt, [], 1, 'omitnan');
    trialmean_galvanic = mean(galvanic_dt, 1, 'omitnan');
    trialstd_galvanic = std(galvanic_dt, [], 1, 'omitnan');
    trialmean_ctrl = mean(ctrl_dt, 1, 'omitnan');
    trialstd_ctrl = std(ctrl_dt, [], 1, 'omitnan');
    trialmean_anodic = mean(anodic_dt, 1, 'omitnan');
    trialstd_anodic = std(anodic_dt, 1, 'omitnan');
    figure;
    set(gca, 'fontsize', 18);
    hold on
    errorbar(pulse_coherences, trialmean_pulse, trialstd_pulse/sqrt(num_trials), 'Color', default_colors(7, :))
    errorbar(galvanic_coherences, trialmean_galvanic, trialstd_galvanic/sqrt(num_trials), 'Color', default_colors(5, :))
    errorbar(control_coherences, trialmean_ctrl, trialstd_ctrl/sqrt(num_trials), 'k')
    errorbar(anodic_coherences, trialmean_anodic, trialstd_anodic/sqrt(num_trials), 'Color', default_colors(6, :))
    scatter(control_coherences, ctrl_dt', 'k', 'filled')
    scatter(pulse_coherences, pulse_dt', [], default_colors(7, :).*ones(length(pulse_dt), 3), 'filled')
    scatter(galvanic_coherences, galvanic_dt', [], default_colors(5, :).*ones(length(galvanic_dt), 3), 'filled')
    scatter(anodic_coherences, anodic_dt', [], default_colors(6, :).*ones(length(anodic_dt), 3), 'filled')
    hold off
    xlabel("Coherence (%)")
    ylabel("Decision Time (s)")
    legend("Pulsatile", "Galvanic", "Control", "Anodic")
    
    %No-decisions
    figure;
    set(gca, 'fontsize', 18);
    hold on
    scatter(pulse_coherences, pulse_nodec*100, [], default_colors(7, :).*ones(length(pulse_nodec), 3), 'filled')
    scatter(galvanic_coherences, galvanic_nodec*100, [], default_colors(5, :).*ones(length(galvanic_nodec), 3), 'filled')
    %scatter(anodic_coherences, anodic_nodec*100, [], default_colors(6, :).*ones(length(anodic_nodec), 3), 'filled')
    scatter(control_coherences, ctrl_nodec*100, 'k', 'filled') 
    hold off
    xlabel("Coherence (%)")
    ylabel("% of trials No-Decision")
end