%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Decision Metrics Accuracy and Decision Time
function plot_decisions(sim_name, pulse_amps, stim_amps, default_colors, num_batch, ...
               pulse_coherences, galvanic_coherences, control_coherences, num_trials)
    %Accuracy
    stim_coeffs = zeros(3, num_batch, 2);
    c = -1:0.01:1;
    stim_ws = zeros(3, num_batch, length(c));
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j <= length(pulse_amps);
        if pulse
            datapath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                [sim_name, stim_amp*1e6]);
            stim_coherences = pulse_coherences;
        elseif stim_amp == 0 %Control
            datapath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e9]);
            stim_coherences = control_coherences;
        else
            datapath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e6]);
            stim_coherences = galvanic_coherences;
        end
        load(strcat(datapath, "/decisions.mat"), "avg_acc", "decisions", "batch_coeffs", "percent_nodec");
        if pulse
            pulse_acc = avg_acc;
            pulse_nodec = percent_nodec;
        elseif stim_amp == 0
            ctrl_acc = avg_acc;
            ctrl_nodec = percent_nodec;
        else
            galvanic_acc = avg_acc;
            galvanic_nodec = percent_nodec;
        end
        [nodec_trial, nodec_c] = find(~decisions);
        for nodec = 1:length(nodec_trial)
            fprintf("c=%0.3f, trial = %0.0f \n", ...
                [stim_coherences(nodec_c(nodec)), nodec_trial(nodec)])
        end
        stim_coeffs(j, :, :) = batch_coeffs;
        for batch = 1:num_batch
            stim_ws(j, batch, :) = logistic_acc(batch_coeffs(batch, :), c);
        end
    end
    pulse_coeffs = reshape(stim_coeffs(1, :, :), [num_batch, 2]);
    galvanic_coeffs = reshape(stim_coeffs(2, :, :), [num_batch, 2]);
    control_coeffs = reshape(stim_coeffs(3, :, :), [num_batch, 2]);
    pulse_w = reshape(stim_ws(1, :, :), [num_batch, length(c)]);
    galvanic_w = reshape(stim_ws(2, :, :), [num_batch, length(c)]);
    control_w = reshape(stim_ws(3, :, :), [num_batch, length(c)]);
    
    figure;
    set(gca, 'fontsize', 18);
    hold on
    scatter(control_coherences, ctrl_acc, 'k', 'filled')
    scatter(pulse_coherences, pulse_acc, [], default_colors(7, :).*ones(length(pulse_acc), 3), 'filled')
    scatter(galvanic_coherences, galvanic_acc, [], default_colors(5, :).*ones(length(galvanic_acc), 3), 'filled')
    plot(c, control_w, "k")
    plot(c, galvanic_w, 'Color', default_colors(5, :))
    plot(c, pulse_w, 'Color', default_colors(7, :))
    hold off
    xlabel("Coherence (%)")
    ylabel("% of trials P1 wins")
    legend("Pulsatile", "Galvanic", "Control")
    
    %Decision Times
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j <= length(pulse_amps);
        if pulse
            datapath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                [sim_name, stim_amp*1e6]);
        elseif stim_amp == 0 %Control
            datapath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e6]);
        else
            datapath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e6]);
        end
        load(strcat(datapath, "/decisions.mat"), "decision_times");
        if pulse
            pulse_dt = decision_times;
        elseif stim_amp == 0
            ctrl_dt = decision_times;
        else
            galvanic_dt = decision_times;
        end
    end
    trialmean_pulse = mean(pulse_dt, 1, 'omitnan');
    trialstd_pulse = std(pulse_dt, [], 1, 'omitnan');
    trialmean_galvanic = mean(galvanic_dt, 1, 'omitnan');
    trialstd_galvanic = std(galvanic_dt, [], 1, 'omitnan');
    trialmean_ctrl = mean(ctrl_dt, 1, 'omitnan');
    trialstd_ctrl = std(ctrl_dt, [], 1, 'omitnan');
    figure;
    set(gca, 'fontsize', 18);
    hold on
    errorbar(pulse_coherences, trialmean_pulse, trialstd_pulse/sqrt(num_trials), 'Color', default_colors(7, :))
    errorbar(galvanic_coherences, trialmean_galvanic, trialstd_galvanic/sqrt(num_trials), 'Color', default_colors(5, :))
    errorbar(control_coherences, trialmean_ctrl, trialstd_ctrl/sqrt(num_trials), 'k')
    scatter(control_coherences, ctrl_dt', 'k', 'filled')
    scatter(pulse_coherences, pulse_dt', [], default_colors(7, :).*ones(length(pulse_dt), 3), 'filled')
    scatter(galvanic_coherences, galvanic_dt', [], default_colors(5, :).*ones(length(pulse_dt), 3), 'filled')
    hold off
    xlabel("Coherence (%)")
    ylabel("Decision Time (s)")
    legend("Pulsatile", "Galvanic", "Control")
    
    %No-decisions
    figure;
    set(gca, 'fontsize', 18);
    hold on
    scatter(pulse_coherences, pulse_nodec*100, [], default_colors(7, :).*ones(length(pulse_nodec), 3), 'filled')
    scatter(galvanic_coherences, galvanic_nodec*100, [], default_colors(5, :).*ones(length(galvanic_nodec), 3), 'filled')
    scatter(control_coherences, ctrl_nodec*100, 'k', 'filled') 
    hold off
    xlabel("Coherence (%)")
    ylabel("% of trials No-Decision")
end