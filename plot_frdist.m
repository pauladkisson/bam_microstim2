%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Firing Rate Distribution over distance from electrode
function plot_frdist(sim_names, ex_c, pulse_amps, stim_amps, t, num_group, num_affected, ...
                     win_start, win_stop, default_colors, ...
                     pulse_coherences, galvanic_coherences, control_coherences, ...
                     anodic_coherences, start_trial, end_trial, num_trials, plot_name)
    for sim_name = sim_names
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
        if contains(sim_name, "Discon")
            discon_frs = stim_frs;
        else
            con_frs = stim_frs;
        end
    end
    pulse_frs = reshape(mean(discon_frs(1, :, :), 2, 'omitnan'), [num_group, 1]);
    galvanic_frs = reshape(mean(discon_frs(2, :, :), 2, 'omitnan'), [num_group, 1]);
    control_frs = reshape(mean(discon_frs(3, :, :), 2, 'omitnan'), [num_group, 1]);
    anodic_frs = reshape(mean(discon_frs(4, :, :), 2, 'omitnan'), [num_group, 1]);
    
    pulse_sems = reshape(std(discon_frs(1, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
    galvanic_sems = reshape(std(discon_frs(2, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
    control_sems = reshape(std(discon_frs(3, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
    anodic_sems = reshape(std(discon_frs(4, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);

    figure;
    set(gca, 'fontsize', 18)
    hold on
    errorbar(ball_rs*1e6, pulse_frs, pulse_sems, '.', ...
        'Color', default_colors(7, :), 'MarkerSize', 20)
    errorbar(ball_rs*1e6, galvanic_frs, galvanic_sems, '.', ...
        'Color', default_colors(5, :), 'MarkerSize', 20)
    errorbar(ball_rs*1e6, control_frs, control_sems, "k.", 'MarkerSize', 20)
    errorbar(ball_rs*1e6, anodic_frs, anodic_sems, '.', ...,
        'MarkerSize', 20, 'Color', default_colors(6, :))
    if contains(plot_name, "zoom")
        ylim([15, 40])
        xlim([0, 2000])
    end
    hold off
    xlabel("Distance from Electrode (um)")
    ylabel("Firing Rate (spk/s)")
    title("Disconnected")
    
    %Full Population Aggregated Activity
    popmean_pulse = reshape(mean(discon_frs(1, :, :), 3, 'omitnan'), [num_trials, 1]);
    popmean_galvanic = reshape(mean(discon_frs(2, :, :), 3, 'omitnan'), [num_trials, 1]);
    popmean_ctrl = reshape(mean(discon_frs(3, :, :), 3, 'omitnan'), [num_trials, 1]);
    popmean_anodic = reshape(mean(discon_frs(4, :, :), 3, 'omitnan'), [num_trials, 1]);
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
    title("Disconnected")
    
    %Connected
    pulse_frs = reshape(mean(con_frs(1, :, :), 2, 'omitnan'), [num_group, 1]);
    galvanic_frs = reshape(mean(con_frs(2, :, :), 2, 'omitnan'), [num_group, 1]);
    control_frs = reshape(mean(con_frs(3, :, :), 2, 'omitnan'), [num_group, 1]);
    anodic_frs = reshape(mean(con_frs(4, :, :), 2, 'omitnan'), [num_group, 1]);
    
    pulse_sems = reshape(std(con_frs(1, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
    galvanic_sems = reshape(std(con_frs(2, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
    control_sems = reshape(std(con_frs(3, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
    anodic_sems = reshape(std(con_frs(4, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);

    figure;
    set(gca, 'fontsize', 18)
    hold on
    errorbar(ball_rs*1e6, pulse_frs, pulse_sems, '.', ...
        'Color', default_colors(7, :), 'MarkerSize', 20)
    errorbar(ball_rs*1e6, galvanic_frs, galvanic_sems, '.', ...
        'Color', default_colors(5, :), 'MarkerSize', 20)
    errorbar(ball_rs*1e6, control_frs, control_sems, "k.", 'MarkerSize', 20)
    errorbar(ball_rs*1e6, anodic_frs, anodic_sems, '.', ...,
        'MarkerSize', 20, 'Color', default_colors(6, :))
    if contains(plot_name, "zoom")
        if any(contains(sim_names, "Con"))
            ylim([15, 40])
            xlim([0, 2000])
        elseif any(contains(sim_names, "Int"))
            ylim([0, 1])
        elseif any(contains(sim_names, "Rec"))
            xlim([0, 2000])
            ylim([40, 70])
        end
    end
    hold off
    xlabel("Distance from Electrode (um)")
    ylabel("Firing Rate (spk/s)")
    title(sim_names(~contains(sim_names, "Discon")))
    
    %Full Population Aggregated Activity
    popmean_pulse = reshape(mean(con_frs(1, :, :), 3, 'omitnan'), [num_trials, 1]);
    popmean_galvanic = reshape(mean(con_frs(2, :, :), 3, 'omitnan'), [num_trials, 1]);
    popmean_ctrl = reshape(mean(con_frs(3, :, :), 3, 'omitnan'), [num_trials, 1]);
    popmean_anodic = reshape(mean(con_frs(4, :, :), 3, 'omitnan'), [num_trials, 1]);
    norm_pulse = popmean_pulse - popmean_ctrl - norm_pulse;
    norm_galvanic = popmean_galvanic - popmean_ctrl - norm_galvanic;
    norm_anodic = popmean_anodic - popmean_ctrl - norm_anodic;
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
    title(sim_names(~contains(sim_names, "Discon")))
    
end