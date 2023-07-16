%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Phase Locking
function plot_phaselock(sim_names, pulse_amps, stim_amps, t, t_task, t_taskoff, t_cut, stim_freq, ...
                        num_group, num_affected, top_N, win_start, win_stop, idx_diff, ...
                        default_colors, start_trial, end_trial, num_trials, ex_c, ...
                        pulse_coherences, galvanic_coherences, control_coherences, anodic_coherences)
    assert(win_start>=t_task & win_start<t_taskoff, "Start of plotting window must be during task")
    assert(win_stop>t_task & win_stop<=t_taskoff, "End of plotting window must be during task")
    pulsetimes = t_task:1/stim_freq:t_taskoff;
    figure(1);
    set(gca, 'Fontsize', 18)
    hold on
    ps_sim_sync = zeros(length(sim_names), num_trials, num_group); 
    for sim_name = sim_names
        stim_sync = zeros(length(stim_amps), num_trials, num_group);
        load(sprintf("Simulation %s/ustim/r.mat", sim_name), "ball_r")
        ball_rs = get_ball_rs(ball_r, num_affected, num_group);
        for j = 1:length(stim_amps)
            c = ex_c(j);
            stim_amp = stim_amps(j);
            pulse = j<=length(pulse_amps);
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
            load(strcat(output_stimpath, "/decisions.mat"), "decisions", "decision_times")
            for trial = start_trial:end_trial
                relative_trial = trial - start_trial + 1;
                if ~contains(sim_name, "Discon") && (...
                        decisions(relative_trial, stim_coherences==c) ~= 1 || ...
                        decision_times(relative_trial, stim_coherences==c) > t_cut)
                    stim_sync(j, relative_trial, :) = NaN;
                    continue %skip trials where P1 doesn't win or decision takes too long
                end
                load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                    "recspikes")
                perc_sync = zeros(num_group, 1);
                for nn = 1:num_group
                    spikeidx = recspikes(int2str(nn));
                    spikeidx = spikeidx(t(spikeidx)>=win_start & t(spikeidx)<win_stop);
                    if isempty(spikeidx)
                        perc_sync(nn) = NaN;
                        continue
                    end
                    for n_spk = 1:length(spikeidx)
                        [t_diff, ~] = min(abs(t(spikeidx(n_spk)) - pulsetimes));
                        if t_diff <= t(idx_diff)
                            perc_sync(nn) = perc_sync(nn) + 1 / length(spikeidx);
                        end
                    end
                end
                stim_sync(j, relative_trial, :) = perc_sync;
            end
        end
        stim_sync = stim_sync * 100;
        pulse_sync = reshape(stim_sync(1, :, :), [num_trials, num_group]);
        galvanic_sync = reshape(stim_sync(2, :, :), [num_trials, num_group]);
        control_sync = reshape(stim_sync(3, :, :), [num_trials, num_group]);
        anodic_sync = reshape(stim_sync(4, :, :), [num_trials, num_group]);
        pulse_trialmean = mean(pulse_sync, 1, 'omitnan');
        galvanic_trialmean = mean(galvanic_sync, 1, 'omitnan');
        control_trialmean = mean(control_sync, 1, 'omitnan');
        anodic_trialmean = mean(anodic_sync, 1, 'omitnan');
        pulse_wins = sum(~isnan(pulse_sync), 1);
        galvanic_wins = sum(~isnan(galvanic_sync), 1);
        control_wins = sum(~isnan(control_sync), 1)
        anodic_wins = sum(~isnan(anodic_sync), 1)
        pulse_sem = std(pulse_sync, [], 1, 'omitnan') ./ sqrt(pulse_wins);
        galvanic_sem = std(galvanic_sync, [], 1, 'omitnan') ./ sqrt(galvanic_wins);
        control_sem =  std(control_sync, [], 1, 'omitnan') ./ sqrt(control_wins);
        anodic_sem = std(anodic_sync, [], 1, 'omitnan') ./ sqrt(anodic_wins);

        if contains(sim_name, "Discon")
            plot_shape = 'o';
        else
            % plot_shape = '^';
            plot_shape = 'o';
        end
        % figure(1);
        figure;
        hold on
        errorbar(ball_rs(1:top_N)*1e6, galvanic_trialmean(1:top_N), galvanic_sem(1:top_N), ...
            plot_shape, 'Color', default_colors(5, :), 'MarkerFaceColor', default_colors(5, :))
        errorbar(ball_rs(1:top_N)*1e6, control_trialmean(1:top_N), control_sem(1:top_N), ...
            plot_shape, 'Color', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
        errorbar(ball_rs(1:top_N)*1e6, pulse_trialmean(1:top_N), pulse_sem(1:top_N), ...
            plot_shape, 'Color', default_colors(7, :), 'MarkerFaceColor', default_colors(7, :))
        errorbar(ball_rs(1:top_N)*1e6, anodic_trialmean(1:top_N), anodic_sem(1:top_N), ...
            plot_shape, 'Color', default_colors(6, :), 'MarkerFaceColor', default_colors(6, :))
        hold off
        xlabel("Distance from Electrode (um)")
        ylabel("Percent of Phaselocked Spikes (%)")
        ylim([0, 100])
        xlim([0, 2000])
        title(sim_name)
        
        % Statistics
        p_ps = ones(num_affected, 1);
        for neuron = 1:num_affected
            ps_sync = reshape(stim_sync(1, :, neuron), [num_trials, 1]);
            ctrl_sync = reshape(stim_sync(3, :, neuron), [num_trials, 1]);
            [~, p_ps(neuron)] = ttest2(ps_sync, ctrl_sync, 'Tail', 'right');
        end
        % Bonferroni correction
        p_ps = p_ps * num_affected;
        sig_thresh = 0.05;
        ps_sig_dist = max(ball_rs(p_ps < sig_thresh));
        fprintf("Simulation %s: PS Synchrony spread to %0.1fum \n", sim_name, ps_sig_dist*1e6)
        ps_sim_sync(sim_name==sim_names, :, :) = pulse_sync;
    end
    % Statistics
    p_sim = ones(num_affected, 1);
    for neuron = 1:num_affected
        con_sync = reshape(ps_sim_sync(1, :, neuron), [num_trials, 1]);
        discon_sync = reshape(ps_sim_sync(2, :, neuron), [num_trials, 1]);
        [~, p_sim(neuron)] = ttest2(con_sync, discon_sync);
    end
    % Bonferroni correction
    p_sim = p_sim * num_affected;
    sig_thresh = 0.05;
    fprintf("Neurons with significantly different synchrony: %0.1f um \n", ...
            ball_rs(p_sim<sig_thresh)*1e6)
end