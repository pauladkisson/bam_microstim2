%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Decision Metrics Accuracy and Decision Time
function plot_decisions(sim_name, pulse_amps, stim_amps, default_colors, brains, ...
                        num_brains, num_batch, ...
                        pulse_coherences, galvanic_coherences, control_coherences)
    %Accuracy
    pulse_acc = zeros(length(brains), length(pulse_coherences));
    galvanic_acc = zeros(length(brains), length(galvanic_coherences));
    ctrl_acc = zeros(1, length(control_coherences));
    pulse_coeffs = zeros(length(brains), 2);
    galvanic_coeffs = zeros(length(brains), 2);
    ctrl_coeffs = zeros(num_batch, 2);
    c = -1:0.01:1;
    w_pulse = zeros(length(brains), length(c));
    w_galvanic = zeros(length(brains), length(c));
    w_ctrl = zeros(num_batch, length(c));
    for brain = brains
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            pulse = j <= length(pulse_amps);
            if pulse
                datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                    [sim_name, brain, stim_amp*1e9]);
                stim_coherences = pulse_coherences;
            elseif stim_amp == 0 %Control
                datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                    [sim_name, brain, stim_amp*1e9]);
                stim_coherences = control_coherences;
                if brain ~= 1
                    continue
                end
            else
                datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                    [sim_name, brain, stim_amp*1e9]);
                stim_coherences = galvanic_coherences;
            end
            load(strcat(datapath, "/decisions.mat"), "avg_final_acc", "final_decisions", "coeffs", "batch_coeffs");
            [nodec_trial, nodec_c] = find(~final_decisions);
            for nodec = 1:length(nodec_trial)
                fprintf("Brain %0.0f, c=%0.3f, trial = %0.0f \n", ...
                    [brain, stim_coherences(nodec_c(nodec)), nodec_trial(nodec)])
            end
            if pulse
                pulse_acc(brain, :) = avg_final_acc;
                pulse_coeffs(brain, :) = coeffs;
                w_pulse(brain, :) = logistic_acc(coeffs, c);
            elseif stim_amp == 0
                ctrl_acc = avg_final_acc;
                ctrl_coeffs = batch_coeffs;
                true_ctrl_coeffs = coeffs;
                for batch = 1:num_batch
                    w_ctrl(batch, :) = logistic_acc(batch_coeffs(batch, :), c);
                end
            else
                galvanic_acc(brain, :) = avg_final_acc;
                galvanic_coeffs(brain, :) = coeffs;
                w_galvanic(brain, :) = logistic_acc(coeffs, c);
            end
        end
    end

    figure;
    hold on
    %errorbar(pulse_coherences, mean(pulse_acc, 1), std(pulse_acc, [], 1)/sqrt(num_brains), 'Color', default_colors(7, :))
    %errorbar(galvanic_coherences, mean(galvanic_acc, 1), std(galvanic_acc, [], 1)/sqrt(num_brains), 'Color', default_colors(5, :))
    %plot(control_coherences, ctrl_acc, 'ko-')
    scatter(control_coherences, ctrl_acc, 'k', 'filled')
    scatter(pulse_coherences, pulse_acc', [], default_colors(7, :).*ones(length(pulse_acc), 3), 'filled')
    scatter(galvanic_coherences, galvanic_acc', [], default_colors(5, :).*ones(length(pulse_acc), 3), 'filled')
    plot(c, mean(w_ctrl, 1), "k")
    plot(c, mean(w_galvanic, 1), 'Color', default_colors(5, :))
    plot(c, mean(w_pulse, 1), 'Color', default_colors(7, :))
    hold off
    xlabel("Coherence (%)")
    ylabel("% of trials P1 wins")
    %legend("Pulsatile", "Galvanic", "Control")

    %Stats
    beta0 = true_ctrl_coeffs(1);
    pulse_beta1 = pulse_coeffs(:, 1) - beta0;
    galvanic_beta1 = galvanic_coeffs(:, 1) - beta0;
    ctrl_beta1 = ctrl_coeffs(:, 1) - beta0;
    pulse_beta_ratio = pulse_beta1 ./ pulse_coeffs(:, 2);
    galvanic_beta_ratio = galvanic_beta1 ./ galvanic_coeffs(:, 2);
    ctrl_beta_ratio = ctrl_beta1 ./ ctrl_coeffs(:, 2);
    groups = [ones(num_brains, 1); 2*ones(num_brains, 1); 3*ones(3, 1)];
    beta_ratios = [pulse_beta_ratio; galvanic_beta_ratio; ctrl_beta_ratio];
    [~, p, stats] = anovan(beta_ratios, groups, 'display', 'off');
    c_val = multcompare(stats, 'display', 'off')
    
    %Decision Times
    pulse_dt = zeros(length(brains), length(pulse_coherences));
    galvanic_dt = zeros(length(brains), length(galvanic_coherences));
    ctrl_dt = zeros(1, length(control_coherences));
    for brain = brains
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            pulse = j <= length(pulse_amps);
            if pulse
                datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                    [sim_name, brain, stim_amp*1e9]);
            elseif stim_amp == 0 %Control
                datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                    [sim_name, brain, stim_amp*1e9]);
                if brain ~= 1
                    continue
                end
            else
                datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                    [sim_name, brain, stim_amp*1e9]);
            end
            load(strcat(datapath, "/decisions.mat"), "final_decision_times", "final_decisions");
            if pulse
                pulse_dt(brain, :) = mean(final_decision_times, 1, 'omitnan');
            elseif stim_amp == 0
                ctrl_dt(brain, :) = mean(final_decision_times, 1, 'omitnan');
            else
                galvanic_dt(brain, :) = mean(final_decision_times, 1, 'omitnan');
            end
        end
    end

    figure;
    hold on
    errorbar(pulse_coherences, mean(pulse_dt, 1), std(pulse_dt, [], 1)/sqrt(num_brains), 'Color', default_colors(7, :))
    errorbar(galvanic_coherences, mean(galvanic_dt, 1), std(galvanic_dt, [], 1)/sqrt(num_brains), 'Color', default_colors(5, :))
    plot(control_coherences, ctrl_dt, 'k')
    scatter(control_coherences, ctrl_dt, 'k', 'filled')
    scatter(pulse_coherences, pulse_dt', [], default_colors(7, :).*ones(length(pulse_dt), 3), 'filled')
    scatter(galvanic_coherences, galvanic_dt', [], default_colors(5, :).*ones(length(pulse_dt), 3), 'filled')
    hold off
    xlabel("Coherence (%)")
    ylabel("Decision Time (s)")
    legend("Pulsatile", "Galvanic", "Control")
    %ylim([0, 1])

    %Stats
    norm_pulse_dt = mean(pulse_dt, 2) - mean(ctrl_dt);
    norm_galvanic_dt = mean(galvanic_dt, 2) - mean(ctrl_dt);
    [~, p] = ttest2(norm_pulse_dt, norm_galvanic_dt)
    norm_pulse_dt_bar = mean(norm_pulse_dt)
    norm_galvanic_dt_bar = mean(norm_galvanic_dt)

    [~, peak_cidx_pulse] = max(pulse_dt, [], 2);
    [~, peak_cidx_galvanic] = max(galvanic_dt, [], 2);
    [~, peak_cidx_ctrl] = max(ctrl_dt, [], 2);
    peak_c_pulse = pulse_coherences(peak_cidx_pulse);
    peak_c_galvanic = galvanic_coherences(peak_cidx_galvanic);
    peak_c_ctrl = control_coherences(peak_cidx_ctrl);
    [~, p] = ttest2(peak_c_pulse, peak_c_galvanic)
    peak_c_pulse_bar = mean(peak_c_pulse)
    peak_c_galvanic_bar = mean(peak_c_galvanic)
    %}
end