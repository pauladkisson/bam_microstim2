%%% Paul Adkisson
%%% 5.15.23
%%% Purpose: Bootstrap decisions to estimate accuracy with confidence

num_bootstraps = 10000;
options = optimoptions('lsqcurvefit', 'Display', 'off');

sim_name = "iScience_Con";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))
default_colors = get(gca, "colororder");

control_coherences = [-100, -51.2, -25.6, -12.8, -6.4, -3.2, 0, 3.2, 6.4, 12.8, 25.6, 51.2, 100] / 100;
pulse_coherences = [-100, -82.6, -69.8, -63.4, -60.2, -57, -53.8, -50.6, ...
                    -44.2, -31.4, -5.8, 0, 43, 100] ./ 100;
galvanic_coherences = [-100, -82.6, -69.8, -63.4, -60.2, -57, -53.8, -50.6, ...
                       -44.2, -31.4, -5.8, 0, 43, 100] ./ 100;
anodic_coherences = fliplr([100, 81.2, 55.6, 42.8, 36.4, 33.2, 30, 26.8, ...
                             23.6, 17.2, 4.4, 0, -21.2, -70, -100]) ./ 100;
pulse_amps = [-10*1e-6];
dc_amps = [-1.4, 0, 1.4]*1e-6;
stim_amps = [pulse_amps, dc_amps];

stim_coeffs = zeros(num_bootstraps, length(stim_amps), 2);
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j <= length(pulse_amps);
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
    load(strcat(output_stimpath, "/decisions.mat"), "decisions", ...
        "decision_times", "avg_acc")
    if pulse
        pulse_acc = avg_acc;
    else
        if stim_amp < 0 %cathodic GS
            galvanic_acc = avg_acc;
        elseif stim_amp == 0 % control
            ctrl_acc = avg_acc;
        else
            anodic_acc = avg_acc;
        end
    end
    num_trials = size(decisions, 1);
    trials = 1:num_trials;
    wait = waitbar(0, "Bootstrapping");
    for bootstrap = 1:num_bootstraps
        waitbar(bootstrap/num_bootstraps, wait, "Bootstrapping");
        trials_resampled = datasample(trials, num_trials);
        decisions_resampled = decisions(trials_resampled, :);
        decision_times_resampled = decision_times(trials_resampled, :);
        avg_acc = zeros(length(stim_coherences), 1);
        for i = 1:length(stim_coherences)
            %omit trials with no decision & those that decide before onset of task-related input
            coherent_decisions = decisions_resampled(decisions_resampled(:, i)~=0 & decision_times_resampled(:, i)>0, i);
            avg_acc(i) = sum(coherent_decisions == 1) / length(coherent_decisions);
        end
        coeffs = lsqcurvefit(@logistic_acc, [1, 1], stim_coherences, avg_acc', [], [], options);
        stim_coeffs(bootstrap, j, :) = coeffs;
    end
end

% Compute Psychometric Parameters
a = reshape(stim_coeffs(:, :, 1), [num_bootstraps, length(stim_amps)]);
b = reshape(stim_coeffs(:, :, 2), [num_bootstraps, length(stim_amps)]);
bias = -a ./ b;
sensitivity = b ./ 4;

figure;
hold on
for j = 1:length(stim_amps)
    histogram(bias(:, j))
end
title("Bias")
legend("Pulsatile", "Galvanic", "Control", "Anodic")

figure;
hold on
for j = 1:length(stim_amps)
    histogram(sensitivity(:, j))
end
title("Sensitivity")
legend("Pulsatile", "Galvanic", "Control", "Anodic")


% Bias Statistics
avg_ps_bias = mean(bias(:, 1));
sem_ps_bias = std(bias(:, 1));
avg_gs_bias = mean(bias(:, 2));
sem_gs_bias = std(bias(:, 2));
avg_ctrl_bias = mean(bias(:, 3));
sem_ctrl_bias = std(bias(:, 3));
avg_an_bias = mean(bias(:, 4));
sem_an_bias = std(bias(:, 4));

bias_cgs_ags = abs(bias(:, 2)) - abs(bias(:, 4));
bias_cgs_ags = bias_cgs_ags(~isnan(bias_cgs_ags));
p_cgs_ags = 2 * sum(bias_cgs_ags < 0) / length(bias_cgs_ags);

bias_ps_cgs = bias(:, 1) - bias(:, 2);
bias_ps_cgs = bias_ps_cgs(~isnan(bias_ps_cgs));
if median(bias_ps_cgs) > 0
    p_ps_cgs = 2 * sum(bias_ps_cgs < 0) / length(bias_ps_cgs);
else
    p_ps_cgs = 2 * sum(bias_ps_cgs > 0) / length(bias_ps_cgs);
end

bias_ctrl = bias(:, 3);
bias_ctrl = bias_ctrl(~isnan(bias_ctrl));
if median(bias_ctrl) > 0
    p_ctrl = 2 * sum(bias_ctrl < 0) / length(bias_ctrl);
else
    p_ctrl = 2 * sum(bias_ctrl > 0) / length(bias_ctrl);
end
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
avg_ps_sens = mean(sensitivity(:, 1));
sem_ps_sens = std(sensitivity(:, 1));
avg_gs_sens = mean(sensitivity(:, 2));
sem_gs_sens = std(sensitivity(:, 2));
avg_ctrl_sens = mean(sensitivity(:, 3));
avg_an_sens = mean(sensitivity(:, 4));
sem_an_sens = std(sensitivity(:, 4));
avg_norm_ps_sens = avg_ps_sens - avg_ctrl_sens;
avg_norm_gs_sens = avg_gs_sens - avg_ctrl_sens;
avg_norm_an_sens = avg_an_sens - avg_ctrl_sens;

sens_ps_ctrl = sensitivity(:, 1) - sensitivity(:, 3);
sens_ps_ctrl = sens_ps_ctrl(~isnan(sens_ps_ctrl));
if median(sens_ps_ctrl) > 0
    p_ps_ctrl = 2 * sum(sens_ps_ctrl < 0) / length(sens_ps_ctrl);
else
    p_ps_ctrl = 2 * sum(sens_ps_ctrl > 0) / length(sens_ps_ctrl);
end

sens_gs_ctrl = sensitivity(:, 2) - sensitivity(:, 3);
sens_gs_ctrl = sens_gs_ctrl(~isnan(sens_gs_ctrl));
if median(sens_gs_ctrl) > 0
    p_gs_ctrl = 2 * sum(sens_gs_ctrl < 0) / length(sens_gs_ctrl);
else
    p_gs_ctrl = 2 * sum(sens_gs_ctrl > 0) / length(sens_gs_ctrl);
end

sens_an_ctrl = sensitivity(:, 4) - sensitivity(:, 3);
sens_an_ctrl = sens_an_ctrl(~isnan(sens_an_ctrl));
if median(sens_an_ctrl) > 0
    p_an_ctrl = 2 * sum(sens_an_ctrl < 0) / length(sens_an_ctrl);
else
    p_an_ctrl = 2 * sum(sens_an_ctrl > 0) / length(sens_an_ctrl);
end
disp("SENSITIVITY")
fprintf([...
    'PS (p=%0.2f) and CGS (p=%0.2f) decreased the sensitivity ', ...
    'by %0.2f +/- %0.2f and %0.2f +/- %0.2f respectively ', ...
    'relative to control. \n'], p_ps_ctrl, p_gs_ctrl, avg_norm_ps_sens, ...
    sem_ps_sens, avg_norm_gs_sens, sem_gs_sens)

fprintf([...
    'AGS (p=%0.2f) increased the sensitivity by %0.2f +/- %0.2f ', ...
    'relative to control. \n'], ...
    p_an_ctrl, avg_norm_an_sens, sem_an_sens)

%%% Figures
c = -1:0.001:1;
bias_median = median(bias, 1);
sensitivity_median = median(sensitivity, 1);
b_median = sensitivity_median .* 4;
a_median = -bias_median .* b_median;
pulse_coeffs = [a_median(1), b_median(1)];
galvanic_coeffs = [a_median(2), b_median(2)];
ctrl_coeffs = [a_median(3), b_median(3)];
anodic_coeffs = [a_median(4), b_median(4)];
bias = sort(bias);
low_idx = floor(size(bias, 1)*0.025);
high_idx = floor(size(bias, 1)*0.975);
bias_low = bias(low_idx, :);
bias_high = bias(high_idx, :);
a_low = -bias_low .* b_median;
a_high = -bias_high .* b_median;
pulse_coeffs_low = [a_low(1), b_median(1)];
pulse_coeffs_high = [a_high(1), b_median(1)];
galvanic_coeffs_low = [a_low(2), b_median(2)];
galvanic_coeffs_high = [a_high(2), b_median(2)];
ctrl_coeffs_low = [a_low(3), b_median(3)];
ctrl_coeffs_high = [a_high(3), b_median(3)];
anodic_coeffs_low = [a_low(4), b_median(4)];
anodic_coeffs_high = [a_high(4), b_median(4)];

c_confidence95 = [c, fliplr(c)];
pulse_w = logistic_acc(pulse_coeffs, c);
pulse_w_low = logistic_acc(pulse_coeffs_low, c);
pulse_w_high = logistic_acc(pulse_coeffs_high, c);
pulse_confidence95 = [pulse_w_low, fliplr(pulse_w_high)];

galvanic_w = logistic_acc(galvanic_coeffs, c);
galvanic_w_low = logistic_acc(galvanic_coeffs_low, c);
galvanic_w_high = logistic_acc(galvanic_coeffs_high, c);
galvanic_confidence95 = [galvanic_w_low, fliplr(galvanic_w_high)];

control_w = logistic_acc(ctrl_coeffs, c);
control_w_low = logistic_acc(ctrl_coeffs_low, c);
control_w_high = logistic_acc(ctrl_coeffs_high, c);
control_confidence95 = [control_w_low, fliplr(control_w_high)];

anodic_w = logistic_acc(anodic_coeffs, c);
anodic_w_low = logistic_acc(anodic_coeffs_low, c);
anodic_w_high = logistic_acc(anodic_coeffs_high, c);
anodic_confidence95 = [anodic_w_low, fliplr(anodic_w_high)];

figure;
set(gca, 'fontsize', 18);
hold on
scatter(control_coherences, ctrl_acc, 100, 'k', 'filled')
scatter(pulse_coherences, pulse_acc, 100, default_colors(7, :), 'filled')
scatter(galvanic_coherences, galvanic_acc, 100, default_colors(5, :), 'filled')
scatter(anodic_coherences, anodic_acc, 100, default_colors(6, :), 'filled')
plot(c, control_w, "k", "Linewidth", 2)
plot(c, galvanic_w, 'Color', default_colors(5, :), "Linewidth", 2)
plot(c, anodic_w, 'Color', default_colors(6, :), "Linewidth", 2)
plot(c, pulse_w, 'Color', default_colors(7, :), "Linewidth", 2)
plot(c, control_w_low, "k")
plot(c, galvanic_w_low, 'Color', default_colors(5, :))
plot(c, anodic_w_low, 'Color', default_colors(6, :))
plot(c, pulse_w_low, 'Color', default_colors(7, :))
plot(c, control_w_high, "k")
plot(c, galvanic_w_high, 'Color', default_colors(5, :))
plot(c, anodic_w_high, 'Color', default_colors(6, :))
plot(c, pulse_w_high, 'Color', default_colors(7, :))
fill(c_confidence95, pulse_confidence95, default_colors(7, :), 'FaceAlpha', 0.25)
fill(c_confidence95, galvanic_confidence95, default_colors(5, :), 'FaceAlpha', 0.25)
fill(c_confidence95, control_confidence95, "k", 'FaceAlpha', 0.25)
fill(c_confidence95, anodic_confidence95, default_colors(6, :), 'FaceAlpha', 0.25)
hold off
xlabel("Coherence (%)")
ylabel("% of trials P1 wins")
legend("Control", "Pulse", "Galvanic", "Anodic")
ylim([0, 1])
title("Bias Confidence Interval")

% Redo for sensitivity confidence interval
sensitivity = sort(sensitivity);
low_idx = floor(size(sensitivity, 1)*0.025);
high_idx = floor(size(sensitivity, 1)*0.975);
sensitivity_low = sensitivity(low_idx, :);
sensitivity_high = sensitivity(high_idx, :);
b_low = sensitivity_low .* 4;
b_high = sensitivity_high .* 4;
a_low = -bias_median .* b_low;
a_high = -bias_median .* b_high;
pulse_coeffs_low = [a_low(1), b_low(1)];
pulse_coeffs_high = [a_high(1), b_high(1)];
galvanic_coeffs_low = [a_low(2), b_low(2)];
galvanic_coeffs_high = [a_high(2), b_high(2)];
ctrl_coeffs_low = [a_low(3), b_low(3)];
ctrl_coeffs_high = [a_high(3), b_high(3)];
anodic_coeffs_low = [a_low(4), b_low(4)];
anodic_coeffs_high = [a_high(4), b_high(4)];

c_confidence95 = [c, fliplr(c)];
pulse_w = logistic_acc(pulse_coeffs, c);
pulse_w_low = logistic_acc(pulse_coeffs_low, c);
pulse_w_high = logistic_acc(pulse_coeffs_high, c);
pulse_confidence95 = [pulse_w_low, fliplr(pulse_w_high)];

galvanic_w = logistic_acc(galvanic_coeffs, c);
galvanic_w_low = logistic_acc(galvanic_coeffs_low, c);
galvanic_w_high = logistic_acc(galvanic_coeffs_high, c);
galvanic_confidence95 = [galvanic_w_low, fliplr(galvanic_w_high)];

control_w = logistic_acc(ctrl_coeffs, c);
control_w_low = logistic_acc(ctrl_coeffs_low, c);
control_w_high = logistic_acc(ctrl_coeffs_high, c);
control_confidence95 = [control_w_low, fliplr(control_w_high)];

anodic_w = logistic_acc(anodic_coeffs, c);
anodic_w_low = logistic_acc(anodic_coeffs_low, c);
anodic_w_high = logistic_acc(anodic_coeffs_high, c);
anodic_confidence95 = [anodic_w_low, fliplr(anodic_w_high)];

figure;
set(gca, 'fontsize', 18);
hold on
scatter(control_coherences, ctrl_acc, 100, 'k', 'filled')
scatter(pulse_coherences, pulse_acc, 100, default_colors(7, :), 'filled')
scatter(galvanic_coherences, galvanic_acc, 100, default_colors(5, :), 'filled')
scatter(anodic_coherences, anodic_acc, 100, default_colors(6, :), 'filled')
plot(c, control_w, "k", "Linewidth", 2)
plot(c, galvanic_w, 'Color', default_colors(5, :), "Linewidth", 2)
plot(c, anodic_w, 'Color', default_colors(6, :), "Linewidth", 2)
plot(c, pulse_w, 'Color', default_colors(7, :), "Linewidth", 2)
plot(c, control_w_low, "k")
plot(c, galvanic_w_low, 'Color', default_colors(5, :))
plot(c, anodic_w_low, 'Color', default_colors(6, :))
plot(c, pulse_w_low, 'Color', default_colors(7, :))
plot(c, control_w_high, "k")
plot(c, galvanic_w_high, 'Color', default_colors(5, :))
plot(c, anodic_w_high, 'Color', default_colors(6, :))
plot(c, pulse_w_high, 'Color', default_colors(7, :))
fill(c_confidence95, pulse_confidence95, default_colors(7, :), 'FaceAlpha', 0.25)
fill(c_confidence95, galvanic_confidence95, default_colors(5, :), 'FaceAlpha', 0.25)
fill(c_confidence95, control_confidence95, "k", 'FaceAlpha', 0.25)
fill(c_confidence95, anodic_confidence95, default_colors(6, :), 'FaceAlpha', 0.25)
hold off
xlabel("Coherence (%)")
ylabel("% of trials P1 wins")
legend("Control", "Pulse", "Galvanic", "Anodic")
ylim([0, 1])
title("Sensitivity Confidence Interval")

% Redo for full parameter space confidence interval
confidence = 0.95;
pulse_psych = logistic_acc_vector(a(:, 1), b(:, 1), c);
pulse_w_low = zeros(1, length(c));
pulse_w_high = zeros(1, length(c));

galvanic_psych = logistic_acc_vector(a(:, 2), b(:, 2), c);
galvanic_w_low = zeros(1, length(c));
galvanic_w_high = zeros(1, length(c));

control_psych = logistic_acc_vector(a(:, 3), b(:, 3), c);
control_w_low = zeros(1, length(c));
control_w_high = zeros(1, length(c));

anodic_psych = logistic_acc_vector(a(:, 4), b(:, 4), c);
anodic_w_low = zeros(1, length(c));
anodic_w_high = zeros(1, length(c));
for i = 1:length(c)
    [low, high] = get_confidence_interval(pulse_psych(:, i), confidence);
    pulse_w_low(i) = low;
    pulse_w_high(i) = high;
    
    [low, high] = get_confidence_interval(galvanic_psych(:, i), confidence);
    galvanic_w_low(i) = low;
    galvanic_w_high(i) = high;
    
    [low, high] = get_confidence_interval(control_psych(:, i), confidence);
    control_w_low(i) = low;
    control_w_high(i) = high;
    
    [low, high] = get_confidence_interval(anodic_psych(:, i), confidence);
    anodic_w_low(i) = low;
    anodic_w_high(i) = high;
end

c_confidence95 = [c, fliplr(c)];
pulse_confidence95 = [pulse_w_low, fliplr(pulse_w_high)];
galvanic_confidence95 = [galvanic_w_low, fliplr(galvanic_w_high)];
control_confidence95 = [control_w_low, fliplr(control_w_high)];
anodic_confidence95 = [anodic_w_low, fliplr(anodic_w_high)];

figure;
set(gca, 'fontsize', 18);
hold on
scatter(control_coherences, ctrl_acc, 100, 'k', 'filled')
scatter(pulse_coherences, pulse_acc, 100, default_colors(7, :), 'filled')
scatter(galvanic_coherences, galvanic_acc, 100, default_colors(5, :), 'filled')
scatter(anodic_coherences, anodic_acc, 100, default_colors(6, :), 'filled')
plot(c, control_w, "k", "Linewidth", 2)
plot(c, galvanic_w, 'Color', default_colors(5, :), "Linewidth", 2)
plot(c, anodic_w, 'Color', default_colors(6, :), "Linewidth", 2)
plot(c, pulse_w, 'Color', default_colors(7, :), "Linewidth", 2)
plot(c, control_w_low, "k")
plot(c, galvanic_w_low, 'Color', default_colors(5, :))
plot(c, anodic_w_low, 'Color', default_colors(6, :))
plot(c, pulse_w_low, 'Color', default_colors(7, :))
plot(c, control_w_high, "k")
plot(c, galvanic_w_high, 'Color', default_colors(5, :))
plot(c, anodic_w_high, 'Color', default_colors(6, :))
plot(c, pulse_w_high, 'Color', default_colors(7, :))
fill(c_confidence95, pulse_confidence95, default_colors(7, :), 'FaceAlpha', 0.25)
fill(c_confidence95, galvanic_confidence95, default_colors(5, :), 'FaceAlpha', 0.25)
fill(c_confidence95, control_confidence95, "k", 'FaceAlpha', 0.25)
fill(c_confidence95, anodic_confidence95, default_colors(6, :), 'FaceAlpha', 0.25)
hold off
xlabel("Coherence (%)")
ylabel("% of trials P1 wins")
legend("Control", "Pulse", "Galvanic", "Anodic")
ylim([0, 1])
title("Full Confidence Interval")

function w = logistic_acc_vector(a, b, c)
    X = a + b.*c;
    w = 1 ./ (1 + exp(-X));
end

function [x_low, x_high] = get_confidence_interval(x, confidence)
    low_percentile = (1 - confidence) / 2;
    high_percentile = 1 - low_percentile;
    x_sorted = sort(x);
    low_idx = floor(size(x, 1)*low_percentile);
    high_idx = floor(size(x, 1)*high_percentile);
    x_low = x_sorted(low_idx);
    x_high = x_sorted(high_idx);
end

