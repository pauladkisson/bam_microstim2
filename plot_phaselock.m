%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Phase Locking
function plot_phaselock(sim_names, pulse_amps, stim_amps, t, t_task, t_taskoff, stim_freq, num_group, ...
                        idx_diff, default_colors, brains, num_brains, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        start_trial, end_trial, num_trials)
    ball_rs = zeros(num_brains, num_group);
    pulsetimes = t_task:1/stim_freq:t_taskoff;
    figure(1);
    hold on
    for sim_name = sim_names
        stim_sync = zeros(length(stim_amps), num_brains, num_group);
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
                load(strcat(output_stimpath, "/decisions.mat"), "final_decisions")
                if ~contains(sim_name, "Discon")
                    num_wins = sum(final_decisions(:, stim_coherences==0)==1, 'all') * ones(num_group, 1);
                else
                    num_wins = num_trials * ones(num_group, 1);
                end
                for c = 0
                    fprintf("coherence: %0.3f%% \n", c*100)
                    nan_neurons = zeros(num_trials, num_group);
                    for trial = start_trial:end_trial
                        relative_trial = trial - start_trial + 1;
                        if ~contains(sim_name, "Discon") && final_decisions(relative_trial, stim_coherences==c) ~= 1
                            continue %skip trials where P1 doesn't win
                        end
                        load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                            "recspikes")
                        perc_sync = zeros(num_group, 1);
                        for nn = 1:num_group
                            spikeidx = recspikes(int2str(nn));
                            spikeidx = spikeidx(t(spikeidx)>=2.5 & t(spikeidx)<t_taskoff);
                            if isempty(spikeidx)
                                perc_sync(nn) = NaN;
                                continue
                            end
                            for n_spk = 1:length(spikeidx)
                                tmp_ts = pulsetimes(t(spikeidx(n_spk))- pulsetimes>=0);
                                [t_diff, ~] = min(abs(t(spikeidx(n_spk)) - tmp_ts));
                                if t_diff <= t(idx_diff)
                                    perc_sync(nn) = perc_sync(nn) + 1 / length(spikeidx);
                                end
                            end
                        end
                        no_spike_neurons = isnan(perc_sync);
                        num_wins(no_spike_neurons) = num_wins(no_spike_neurons) - 1;
                        perc_sync(no_spike_neurons) = [];
                        stim_sync(j, brain, ~no_spike_neurons) = ...
                            reshape(stim_sync(j, brain, ~no_spike_neurons), size(perc_sync)) + perc_sync;
                        nan_neurons(trial, :) = no_spike_neurons;
                    end
                    if any(all(nan_neurons==1, 1)) %no spikes for all 36 trials
                        stim_sync(j, brain, all(nan_neurons==1, 1)) = NaN;
                    end
                end
                stim_sync(j, brain, :) = reshape(stim_sync(j, brain, :), size(num_wins)) ./ num_wins;
            end
        end
        nan_sync = isnan(stim_sync); 
        stim_sync = stim_sync * 100;
        pulse_sync = reshape(stim_sync(1, :, :), [num_brains, num_group]);
        pulse_nan = reshape(nan_sync(1, :, :), [num_brains, num_group]);
        galvanic_sync = reshape(stim_sync(2, :, :), [num_brains, num_group]);
        galvanic_nan = reshape(nan_sync(2, :, :), [num_brains, num_group]);
        control_sync = reshape(stim_sync(3, 1, :), [1, num_group]);
        control_nan = reshape(nan_sync(3, 1, :), [1, num_group]);

        top_N = num_group;
        total_N = top_N*num_brains;
        tot_pulse_sync = reshape(pulse_sync(:, 1:top_N), [total_N, 1]);
        tot_galvanic_sync = reshape(galvanic_sync(:, 1:top_N), [total_N, 1]);
        tot_ctrl_sync = reshape(control_sync(1:top_N), [top_N, 1]);
        tot_ball_rs = reshape(ball_rs(:, 1:top_N), [total_N, 1]);
        tot_ctrl_ball_r = ctrl_ball_r(1:top_N);
        tot_pulse_nan = reshape(pulse_nan(:, 1:top_N), [total_N, 1]);
        tot_galvanic_nan = reshape(galvanic_nan(:, 1:top_N), [total_N, 1]);
        tot_ctrl_nan = reshape(control_nan(1:top_N), [top_N, 1]);

        figure(1);
        %figure;
        hold on
        if contains(sim_name, "Discon")
            scatter(tot_ball_rs(~tot_pulse_nan)*1e6, tot_pulse_sync(~tot_pulse_nan), [], default_colors(7, :), 'filled')
            scatter(tot_ball_rs(~tot_galvanic_nan)*1e6, tot_galvanic_sync(~tot_galvanic_nan), [], default_colors(5, :), 'filled')
            scatter(tot_ctrl_ball_r(~tot_ctrl_nan)*1e6, tot_ctrl_sync(~tot_ctrl_nan), [], "k", 'filled')
        else
            scatter(tot_ball_rs(~tot_pulse_nan)*1e6, tot_pulse_sync(~tot_pulse_nan), [], default_colors(7, :), 'filled', '^')
            scatter(tot_ball_rs(~tot_galvanic_nan)*1e6, tot_galvanic_sync(~tot_galvanic_nan), [], default_colors(5, :), 'filled', '^')
            scatter(tot_ctrl_ball_r(~tot_ctrl_nan)*1e6, tot_ctrl_sync(~tot_ctrl_nan), [], "k", 'filled', '^')
        end
        xlabel("Distance from Electrode (um)")
        ylabel("Synchronized APs to Pulse Times (%)")

        %trial/pop means
        pulse_trialpop_mean = mean(pulse_sync, 2, 'omitnan');
        galvanic_trialpop_mean = mean(galvanic_sync, 2, 'omitnan');
        %ctrl_mean = mean(control_sync, 'omitnan');
        %trl_std = std(control_sync, 'omitnan');
        ctrl_mean = 11.3484; %use connected values
        ctrl_std = 1.8900;
        pulse_norm_mean = pulse_trialpop_mean - ctrl_mean;
        galvanic_norm_mean = galvanic_trialpop_mean - ctrl_mean;
        stim_means = [mean(galvanic_norm_mean); mean(pulse_norm_mean)];
        stim_sems = [std(galvanic_norm_mean)/sqrt(num_brains); std(pulse_norm_mean)/sqrt(num_brains)];

        [~, p] = ttest2(galvanic_norm_mean, pulse_norm_mean)
        galvanic_norm_bar = mean(galvanic_norm_mean)
        galvanic_norm_sem = std(galvanic_norm_mean) ./ sqrt(length(galvanic_norm_mean))
        pulse_norm_bar = mean(pulse_norm_mean)
        pulse_norm_sem = std(pulse_norm_mean) ./ sqrt(length(pulse_norm_mean))

        if ~contains(sim_name, "Discon")
            connected_stim_means = stim_means;
            connected_gs_norm_mean = galvanic_norm_mean;
            connected_ps_norm_mean = pulse_norm_mean;
            connected_ps_tot = tot_pulse_sync;
        else
            disconnected_stim_means = stim_means;
            disconnected_gs_norm_mean = galvanic_norm_mean;
            disconnected_ps_norm_mean = pulse_norm_mean;
            disconnected_ps_tot = tot_pulse_sync;
        end
        % %affected
        pulse_brainsync = reshape(stim_sync(1, :, :), [num_brains, num_group]);
        galvanic_brainsync = reshape(stim_sync(2, :, :), [num_brains, num_group]);
        pulse_affected = pulse_brainsync > ctrl_mean + ctrl_std*3 | pulse_brainsync < ctrl_mean - ctrl_std*3;
        disp(sim_name)
        ctrl_mean
        ctrl_std
        affected_distance = max(ball_rs(pulse_affected), [], 'all')*1e6
        galvanic_affected = galvanic_brainsync > ctrl_mean + ctrl_std*3 | galvanic_brainsync < ctrl_mean - ctrl_std*3;
        pulse_percent_affected = sum(pulse_affected, 2) / num_group;
        galvanic_percent_affected = sum(galvanic_affected, 2) / num_group;
        avg_pulse_percent_affected = mean(pulse_percent_affected)
        sem_pulse_percent_affected = std(pulse_percent_affected) / sqrt(length(pulse_percent_affected))
    end
    legend(["Disconnected", "", "", "Connected"])

    hold off

    barcolor = [default_colors(5, :);
                default_colors(7, :)];
    dis_pos = [1, 2];
    con_pos = [4, 5];
    figure;
    hold on
    b = bar(dis_pos, disconnected_stim_means, ...
        'FaceColor', 'flat', 'CData', barcolor);
    b = bar(con_pos, connected_stim_means, ...
        'FaceColor', 'flat', 'CData', barcolor);
    ylabel("Change in Synchronized APs to Pulse Times (%)")
    plot([ones(1, num_brains)*dis_pos(1); ones(1, num_brains)*dis_pos(2)], [disconnected_gs_norm_mean'; disconnected_ps_norm_mean'], 'ko-')
    plot([ones(1, num_brains)*con_pos(1); ones(1, num_brains)*con_pos(2)], [connected_gs_norm_mean'; connected_ps_norm_mean'], 'ko-')
    hold off
    xticks([mean(dis_pos), mean(con_pos)])
    xticklabels(["Dis", "Con"])
    ylim([0, 15])

    [~, ~, stats] = anova1([connected_gs_norm_mean, disconnected_gs_norm_mean, connected_ps_norm_mean, disconnected_ps_norm_mean]);
    c = multcompare(stats)

    [~, p] = kstest2(connected_ps_tot, disconnected_ps_tot)
end