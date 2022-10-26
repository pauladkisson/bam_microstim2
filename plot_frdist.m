%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Firing Rate Distribution over distance from electrode
function plot_frdist(sim_name, ex_c, pulse_amps, stim_amps, t, num_group, num_affected, ...
                     win_start, win_stop, default_colors, ...
                     pulse_coherences, galvanic_coherences, control_coherences, ...
                     start_trial, end_trial, num_trials, plot_name)
    stim_frs = zeros(length(stim_amps), num_trials, num_group);
    load(sprintf("Simulation %s/ustim/r.mat", sim_name), "ball_r")
    ball_rs = get_ball_rs(ball_r, num_affected, num_group);
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j<=length(pulse_amps);
        if pulse
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                [sim_name, stim_amp*1e6]);
            stim_coherences = pulse_coherences;
        else
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e6]);
            if stim_amp == 0
                stim_coherences = control_coherences;
            else
                stim_coherences = galvanic_coherences;
            end
        end
        if sim_name == "EMBC Disconnected"
            stim_coherences = 0;
        end
        load(strcat(output_stimpath, "/decisions.mat"), "decisions")
        for trial = start_trial:end_trial
            relative_trial = trial - start_trial + 1;
            if (plot_name == "p1_wins" && decisions(relative_trial, stim_coherences==ex_c) ~= 1) || ...
                    (plot_name == "p1_loses" && decisions(relative_trial, stim_coherences==ex_c) ~= 2)
                stim_frs(j, relative_trial, :) = NaN;
                continue
            end
            load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, trial])), ...
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
    ctrl_mean = mean(control_frs);
    ctrl_std = std(control_frs);

    figure;
    set(gca, 'fontsize', 18)
    hold on
    scatter(ball_rs*1e6, pulse_frs, [], ones(num_group, 3).*default_colors(7, :), 'filled')
    scatter(ball_rs*1e6, galvanic_frs, [], ones(num_group, 3).*default_colors(5, :), 'filled')
    scatter(ball_rs*1e6, control_frs, [], "k", 'filled')
    yline(ctrl_mean+3*ctrl_std, "k--")
    yline(ctrl_mean-3*ctrl_std, "k--") 
    hold off
    xlabel("Distance from Electrode (um)")
    ylabel("Firing Rate (spk/s)")
    if contains(sim_name, "Discon")
        title("Disconnected")
    else
        title("Connected")
    end
    
    %Full Population Aggregated Activity
    popmean_pulse = reshape(mean(stim_frs(1, :, :), 3, 'omitnan'), [num_trials, 1]);
    popmean_galvanic = reshape(mean(stim_frs(2, :, :), 3, 'omitnan'), [num_trials, 1]);
    norm_pulse = popmean_pulse - ctrl_mean;
    norm_galvanic = popmean_galvanic - ctrl_mean;
    stim_means = [mean(norm_galvanic), mean(norm_pulse)];
    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(stim_means);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(7, :)];
    plot([ones(1, num_trials); 2*ones(1, num_trials)], [norm_galvanic'; norm_pulse'], 'ko-')
    hold off
    xticks([1, 2])
    xticklabels(["Galvanic", "Pulsatile"])
    ylabel("Change in Firing Rate (spk/s)")
    if contains(sim_name, "Discon")
        title("Disconnected")
    else
        title("Connected")
    end
    
    %Affected P1 Aggregated
    popmean_pulse = reshape(mean(stim_frs(1, :, 1:num_affected), 3, 'omitnan'), [num_trials, 1]);
    popmean_galvanic = reshape(mean(stim_frs(2, :, 1:num_affected), 3, 'omitnan'), [num_trials, 1]);
    norm_pulse = popmean_pulse - ctrl_mean;
    norm_galvanic = popmean_galvanic - ctrl_mean;
    stim_means = [mean(norm_galvanic), mean(norm_pulse)];
    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(stim_means);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(7, :)];
    plot([ones(1, num_trials); 2*ones(1, num_trials)], [norm_galvanic'; norm_pulse'], 'ko-')
    hold off
    xticks([1, 2])
    xticklabels(["Galvanic", "Pulsatile"])
    ylabel("Change in Firing Rate (spk/s)")
    if contains(sim_name, "Discon")
        title("Disconnected")
    else
        title("Connected")
    end
    
    %Unaffected P1 Aggregated
    popmean_pulse = reshape(mean(stim_frs(1, :, num_affected+1:end), 3, 'omitnan'), [num_trials, 1]);
    popmean_galvanic = reshape(mean(stim_frs(2, :, num_affected+1:end), 3, 'omitnan'), [num_trials, 1]);
    norm_pulse = popmean_pulse - ctrl_mean;
    norm_galvanic = popmean_galvanic - ctrl_mean;
    stim_means = [mean(norm_galvanic), mean(norm_pulse)];
    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(stim_means);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(7, :)];
    plot([ones(1, num_trials); 2*ones(1, num_trials)], [norm_galvanic'; norm_pulse'], 'ko-')
    hold off
    xticks([1, 2])
    xticklabels(["Galvanic", "Pulsatile"])
    ylabel("Change in Firing Rate (spk/s)")
    if contains(sim_name, "Discon")
        title("Disconnected")
    else
        title("Connected")
    end
end