%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Firing Rate Distribution over distance from electrode
function plot_frdist(sim_name, ex_c, pulse_amps, stim_amps, t, num_group, win_start, ...
                     win_stop, default_colors, brains, num_brains, ...
                     pulse_coherences, galvanic_coherences, control_coherences, ...
                     start_trial, end_trial, num_trials, plot_name) 
    stim_frs = zeros(length(stim_amps), num_brains, num_group);
    ball_rs = zeros(num_brains, num_group);
    for brain = brains
        fprintf("brain %0.0f \n", brain)
        load(sprintf("Simulation %s/brain%0.0f/r.mat", [sim_name, brain]), "ball_r")
        ball_rs(brain, :) = ball_r;
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            pulse = j<=length(pulse_amps);
            if pulse
                output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                    [sim_name, brain, stim_amp*1e9]);
                stim_coherences = pulse_coherences;
            else
                output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                    [sim_name, brain, stim_amp*1e9]);
                if stim_amp == 0
                    stim_coherences = control_coherences;
                    ctrl_ball_r = ball_r;
                    if brain ~= 1
                        continue
                    end
                else
                    stim_coherences = galvanic_coherences;
                end
            end
            if sim_name == "EMBC Disconnected"
                stim_coherences = 0;
            end
            load(strcat(output_stimpath, "/decisions.mat"), "final_decisions")
            if plot_name == "p1_wins"
                num_wins = sum(final_decisions==1, 'all');
                analyze_coherences = stim_coherences;
            elseif plot_name == "p1_loses"
                num_wins = sum(final_decisions==2, 'all');
                analyze_coherences = stim_coherences;
            else
                num_wins = num_trials;
                analyze_coherences = ex_c;
            end
            for ex_c = analyze_coherences
                for trial = start_trial:end_trial
                    relative_trial = trial - start_trial + 1;
                    if (plot_name == "p1_wins" && final_decisions(relative_trial, stim_coherences==ex_c) ~= 1) || ...
                            (plot_name == "p1_loses" && final_decisions(relative_trial, stim_coherences==ex_c) ~= 2)
                        continue
                    end
                    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, trial])), ...
                        "recspikes")
                    g1_taskfrs = zeros(num_group, 1);
                    for nn = 1:num_group
                        spiketimes = t(recspikes(int2str(nn)));
                        g1_taskfrs(nn) = sum(spiketimes>=win_start & spiketimes<win_stop) / (win_stop - win_start);
                    end
                    stim_frs(j, brain, :) = reshape(stim_frs(j, brain, :), size(g1_taskfrs)) + g1_taskfrs ./ num_wins;
                end
            end
        end
    end
    ctrl_mean = mean(stim_frs(3, 1, :));
    ctrl_std = std(stim_frs(3, 1, :));
    %top_N = num_group * 0.2; %show top 20% closest to the electrode
    top_N = num_group;
    total_N = top_N*num_brains;
    pulse_frs = reshape(stim_frs(1, :, 1:top_N), [total_N, 1]);
    galvanic_frs = reshape(stim_frs(2, :, 1:top_N), [total_N, 1]);
    ctrl_frs = reshape(stim_frs(3, 1, 1:top_N), [top_N, 1]);
    ball_rs = reshape(ball_rs(:, 1:top_N), [total_N, 1]);

    figure;
    set(gca, 'fontsize', 18)
    hold on
    scatter(ball_rs*1e6, pulse_frs, [], ones(total_N, 3).*default_colors(7, :), 'filled')
    scatter(ball_rs*1e6, galvanic_frs, [], ones(total_N, 3).*default_colors(5, :), 'filled')
    scatter(ctrl_ball_r(1:top_N)*1e6, ctrl_frs, [], "k", 'filled')
    yline(ctrl_mean+3*ctrl_std, "k--")
    yline(ctrl_mean-3*ctrl_std, "k--") 
    hold off
    xlabel("Distance from Electrode (um)")
    ylabel("Firing Rate (spk/s)")
    if sim_name == "EMBC Disconnected" || sim_name == "DepolBlockDiscon" || sim_name == "Test"
        title("Disconnected")
    else
        title("Connected")
    end

    %trial/pop means
    pulse_trialpop_mean = mean(reshape(stim_frs(1, :, :), [num_brains, num_group]), 2);
    galvanic_trialpop_mean = mean(reshape(stim_frs(2, :, :), [num_brains, num_group]), 2);
    pulse_norm_mean = pulse_trialpop_mean - ctrl_mean;
    galvanic_norm_mean = galvanic_trialpop_mean - ctrl_mean;
    stim_means = [mean(galvanic_norm_mean); mean(pulse_norm_mean)];
    stim_sems = [std(galvanic_norm_mean)/sqrt(num_brains); std(pulse_norm_mean)/sqrt(num_brains)];

    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(stim_means);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(7, :)];
    plot([ones(1, num_brains); 2*ones(1, num_brains)], [galvanic_norm_mean'; pulse_norm_mean'], 'ko-')
    hold off
    xticks([1, 2])
    xticklabels(["Galvanic", "Pulsatile"])
    ylabel("Change in Firing Rate (spk/s)")
    if sim_name == "EMBC Disconnected" || sim_name == "DepolBlockDiscon" || sim_name == "Test"
        title("Disconnected")
    else
        title("Connected")
    end

    % %affected Note: only considering increasing firing rate as affected
    pulse_brainfrs = reshape(stim_frs(1, :, :), [num_brains, num_group]);
    galvanic_brainfrs = reshape(stim_frs(2, :, :), [num_brains, num_group]);
    pulse_affected = pulse_brainfrs > ctrl_mean + ctrl_std*3;
    galvanic_affected = galvanic_brainfrs > ctrl_mean + ctrl_std*3;
    pulse_percent_affected = sum(pulse_affected, 2) / size(pulse_affected, 2)
    galvanic_percent_affected = sum(galvanic_affected, 2) / size(galvanic_affected, 2)
    [~, p_perc_act] = ttest2(pulse_percent_affected, galvanic_percent_affected)

    % %of neurons affected by pulse by also affected by galvanic
    num_affected = 0;
    for brain = brains
        for nn = 1:num_group
            if pulse_affected(brain, nn) && galvanic_affected(brain, nn)
                num_affected = num_affected + 1;
            end
        end
    end
    percent_pulse_affected_affected_by_galvanic = num_affected / sum(pulse_affected, 'all');

    %stats
    [~, p_ttest] = ttest2(galvanic_norm_mean, pulse_norm_mean)
    galvanic_norm_bar = mean(galvanic_norm_mean)
    pulse_norm_bar = mean(pulse_norm_mean)
    galvanic_norm_sem = std(galvanic_norm_mean) ./ sqrt(length(galvanic_norm_mean))
    pulse_norm_sem = std(pulse_norm_mean) ./ sqrt(length(pulse_norm_mean))
    pulse_trialbrainmean = mean(reshape(stim_frs(1, :, :), [num_brains, num_group]), 1);
    galvanic_trialbrainmean = mean(reshape(stim_frs(2, :, :), [num_brains, num_group]), 1);
    figure;
    hold on
    [~, p_kstest] = kstest2(pulse_frs, galvanic_frs)
end