%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Firing Rate Distribution over distance from electrode
function plot_frdist(sim_name, ex_c, pulse_amps, stim_amps, t, num_group, num_affected, ...
                     win_start, win_stop, default_colors, ...
                     pulse_coherences, galvanic_coherences, control_coherences, ...
                     anodic_coherences, start_trial, end_trial, num_trials, plot_name)
    stim_frs = zeros(length(stim_amps), num_trials, num_group);
    load(sprintf("Simulation %s/ustim/r.mat", sim_name), "ball_r")
    ball_rs = get_ball_rs(ball_r, num_affected, num_group);
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j<=length(pulse_amps);
        c = ex_c(j);
        if pulse
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                [sim_name, stim_amp*1e6]);
            stim_coherences = pulse_coherences;
        else
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e6]);
            if stim_amp < 0 %cathodic GS
               stim_coherences = galvanic_coherences;
            elseif stim_amp == 0
                stim_coherences = control_coherences;
            else %anodic GS
                stim_coherences = anodic_coherences;
            end
        end
        load(strcat(output_stimpath, "/decisions.mat"), "decisions")
        for trial = start_trial:end_trial
            relative_trial = trial - start_trial + 1;
            if (plot_name == "p1_wins" && decisions(relative_trial, stim_coherences==c) ~= 1) || ...
                    (plot_name == "p1_loses" && decisions(relative_trial, stim_coherences==c) ~= 2)
                stim_frs(j, relative_trial, :) = NaN;
                continue
            end
            load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                "recspikes")
            g1_taskfrs = zeros(num_group, 1);
            for nn = 1:num_group
                spiketimes = t(recspikes(int2str(nn)));
                g1_taskfrs(nn) = sum(spiketimes>=win_start & spiketimes<win_stop) / (win_stop - win_start);
            end
            stim_frs(j, trial, :) = g1_taskfrs;
        end
    end
    pulse_frs = reshape(mean(stim_frs(1, :, :), 2, 'omitnan'), [num_group, 1]);
    galvanic_frs = reshape(mean(stim_frs(2, :, :), 2, 'omitnan'), [num_group, 1]);
    control_frs = reshape(mean(stim_frs(3, :, :), 2, 'omitnan'), [num_group, 1]);
    %anodic_frs = reshape(mean(stim_frs(4, :, :), 2, 'omitnan'), [num_group, 1]);

    figure;
    set(gca, 'fontsize', 18)
    hold on
    scatter(ball_rs*1e6, pulse_frs, [], ones(num_group, 3).*default_colors(7, :), 'filled')
    scatter(ball_rs*1e6, galvanic_frs, [], ones(num_group, 3).*default_colors(5, :), 'filled')
    scatter(ball_rs*1e6, control_frs, [], "k", 'filled')
    %scatter(ball_rs*1e6, anodic_frs, [], ones(num_group, 3).*default_colors(6, :), 'filled')
    ctrl_mean = mean(stim_frs(3, :, :), 'all');
    ctrl_std = std(mean(stim_frs(3, :, :), 2), [], 'all');
    yline(ctrl_mean+3*ctrl_std, 'k--')
    hold off
    xlabel("Distance from Electrode (um)")
    ylabel("Firing Rate (spk/s)")
    if contains(sim_name, "Discon")
        title("Disconnected")
    else
        title("Connected")
    end
    %percent activated
    ctrl_mean = mean(stim_frs(3, :, :), 'all');
    ctrl_std = std(mean(stim_frs(3, :, :), 2), [], 'all');
    ps_perc_acts = sum(stim_frs(1, :, :) > (ctrl_mean + 3*ctrl_std), 3) / num_affected;
    gs_perc_acts = sum(stim_frs(2, :, :) > (ctrl_mean + 3*ctrl_std), 3) / num_affected;
    %ps_perc_acts = sum(pulse_frs > (ctrl_mean + 3*ctrl_std)) / num_affected;
    %gs_perc_acts = sum(galvanic_frs > (ctrl_mean + 3*ctrl_std)) / num_affected;
    mean_ps_act = mean(ps_perc_acts);
    sem_ps_act = std(ps_perc_acts) / sqrt(num_trials);
    mean_gs_act = mean(gs_perc_acts);
    sem_gs_act = std(gs_perc_acts) / sqrt(num_trials);
    fprintf("PS affects %0.2f%% +/-%0.2f%% of P1 neurons \n", [mean_ps_act*100, sem_ps_act*100]);
    fprintf("GS affects %0.2f%% +/-%0.2f%% of P1 neurons \n", [mean_gs_act*100, sem_gs_act*100]);
    
    %Errorbars in space
    pulse_stds = reshape(std(stim_frs(1, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
    galvanic_stds = reshape(std(stim_frs(2, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
    control_stds = reshape(std(stim_frs(3, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
    figure;
    set(gca, 'fontsize', 18)
    hold on
    errorbar(ball_rs*1e6, pulse_frs, pulse_stds, 'Color', default_colors(7, :))
    errorbar(ball_rs*1e6, galvanic_frs, galvanic_stds, 'Color', default_colors(5, :))
    errorbar(ball_rs*1e6, control_frs, control_stds, "k")
    %scatter(ball_rs*1e6, anodic_frs, [], ones(num_group, 3).*default_colors(6, :), 'filled')
    ctrl_mean = mean(stim_frs(3, :, :), 'all');
    ctrl_std = std(mean(stim_frs(3, :, :), 2), [], 'all');
    yline(ctrl_mean+3*ctrl_std, 'k--')
    hold off
    xlabel("Distance from Electrode (um)")
    ylabel("Firing Rate (spk/s)")
    if contains(sim_name, "Discon")
        title("Disconnected")
    else
        title("Connected")
    end
    %{
    %Full Population Aggregated Activity
    popmean_pulse = reshape(mean(stim_frs(1, :, :), 3, 'omitnan'), [num_trials, 1]);
    popmean_galvanic = reshape(mean(stim_frs(2, :, :), 3, 'omitnan'), [num_trials, 1]);
    popmean_ctrl = reshape(mean(stim_frs(3, :, :), 3, 'omitnan'), [num_trials, 1]);
    popmean_anodic = reshape(mean(stim_frs(4, :, :), 3, 'omitnan'), [num_trials, 1]);
    norm_pulse = popmean_pulse - popmean_ctrl;
    norm_galvanic = popmean_galvanic - popmean_ctrl;
    norm_anodic = popmean_anodic - popmean_ctrl;
    stim_means = [mean(norm_galvanic, 'omitnan'), mean(norm_anodic, 'omitnan'), ...
                  mean(norm_pulse, 'omitnan')];
    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(stim_means);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(6, :); default_colors(7, :)];
    x = [ones(1, num_trials); 2*ones(1, num_trials); 3*ones(1, num_trials)];
    y = [norm_galvanic'; norm_anodic'; norm_pulse'];
    plot(x, y, 'ko')
    hold off
    xticks([1, 2, 3])
    xticklabels(["Galvanic", "Anodic", "Pulsatile"])
    ylabel("Change in Firing Rate (spk/s)")
    title("Full P1")
    
    %Affected P1 Aggregated
    popmean_pulse = reshape(mean(stim_frs(1, :, 1:num_affected), 3, 'omitnan'), [num_trials, 1]);
    popmean_galvanic = reshape(mean(stim_frs(2, :, 1:num_affected), 3, 'omitnan'), [num_trials, 1]);
    popmean_ctrl = reshape(mean(stim_frs(3, :, 1:num_affected), 3, 'omitnan'), [num_trials, 1]);
    popmean_anodic = reshape(mean(stim_frs(4, :, 1:num_affected), 3, 'omitnan'), [num_trials, 1]);
    norm_pulse = popmean_pulse - popmean_ctrl;
    norm_galvanic = popmean_galvanic - popmean_ctrl;
    norm_anodic = popmean_anodic - popmean_ctrl;
    stim_means = [mean(norm_galvanic, 'omitnan'), mean(norm_anodic, 'omitnan'), ...
                  mean(norm_pulse, 'omitnan')];
    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(stim_means);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(6, :); default_colors(7, :)];
    x = [ones(1, num_trials); 2*ones(1, num_trials); 3*ones(1, num_trials)];
    y = [norm_galvanic'; norm_anodic'; norm_pulse'];
    plot(x, y, 'ko')
    hold off
    xticks([1, 2, 3])
    xticklabels(["Galvanic", "Anodic", "Pulsatile"])
    ylabel("Change in Firing Rate (spk/s)")
    title("P1 Affected")
    
    %Unaffected P1 Aggregated
    popmean_pulse = reshape(mean(stim_frs(1, :, num_affected+1:end), 3, 'omitnan'), [num_trials, 1]);
    popmean_galvanic = reshape(mean(stim_frs(2, :, num_affected+1:end), 3, 'omitnan'), [num_trials, 1]);
    popmean_ctrl = reshape(mean(stim_frs(3, :, num_affected+1:end), 3, 'omitnan'), [num_trials, 1]);
    popmean_anodic = reshape(mean(stim_frs(4, :, num_affected+1:end), 3, 'omitnan'), [num_trials, 1]);
    norm_pulse = popmean_pulse -popmean_ctrl;
    norm_galvanic = popmean_galvanic - popmean_ctrl;
    norm_anodic = popmean_anodic - popmean_ctrl;
    stim_means = [mean(norm_galvanic, 'omitnan'), mean(norm_anodic, 'omitnan'), ...
              mean(norm_pulse, 'omitnan')];
    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(stim_means);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(6, :); default_colors(7, :)];
    x = [ones(1, num_trials); 2*ones(1, num_trials); 3*ones(1, num_trials)];
    y = [norm_galvanic'; norm_anodic'; norm_pulse'];
    plot(x, y, 'ko')
    hold off
    xticks([1, 2, 3])
    xticklabels(["Galvanic", "Anodic", "Pulsatile"])
    ylabel("Change in Firing Rate (spk/s)")
    title("P1 Unaffected")
    %}
end