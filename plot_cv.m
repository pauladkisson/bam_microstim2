%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Coefficient of Variation (CV)
function plot_cv(sim_name, sim_names, pulse_amps, stim_amps, t, t_cut, N, top_N, num_group, ...
                 num_affected, win_size, cv_window, default_colors, ex_c, ex_trial, ...
                 ex_neuron, pulse_coherences, galvanic_coherences, control_coherences, ...
                 anodic_coherences, start_trial, end_trial, num_trials, plot_name)
    if plot_name == "ex_neuron" || plot_name == "ex_trial"
        figure;
        set(gca, 'Fontsize', 18);
        hold on
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            pulse = j<=length(pulse_amps);
            if pulse
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                    [sim_name, stim_amp*1e6]);
                stim_color = default_colors(7, :);
            else
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                    [sim_name, stim_amp*1e6]);
                if stim_amp < 0 %cathodic GS
                    stim_color = default_colors(5, :);
                elseif stim_amp == 0 %Control
                    stim_color = [0, 0, 0];
                else %Anodic GS
                    stim_color = default_colors(6, :);
                end
            end
            load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), ...
                "recspikes")
            isi = zeros(length(t), N);
            for nn = 1:N
                spikeidx = recspikes(int2str(nn));
                if length(spikeidx) < 2
                    isi(:, nn) = NaN;
                    continue
                end
                neuron_isi_idx = diff(spikeidx);
                pre_isi_idx = max(neuron_isi_idx(1), spikeidx(1));
                isi(1:spikeidx(1)-1, nn) = t(pre_isi_idx);
                for i = 1:length(neuron_isi_idx)-1
                    isi_idx = neuron_isi_idx(i);
                    spikei = spikeidx(i);
                    next_spikei = spikeidx(i+1);
                    isi(spikei:next_spikei-1, nn) = t(isi_idx);
                end
                post_isi_idx = max(neuron_isi_idx(end), length(t)-spikeidx(end));
                isi(spikeidx(end):end, nn) = t(post_isi_idx);
            end
            p1_isi = isi(cv_window, 1:num_group);
            cv = movstd(p1_isi, win_size, [], 1) ./ movmean(p1_isi, win_size, 1);
            omit_idx = find(any(isnan(cv), 1));
            new_ex_neuron = ex_neuron;
            for idx = omit_idx
                if ex_neuron > omit_idx
                    new_ex_neuron = new_ex_neuron-1;
                elseif ex_neuron == omit_idx
                    disp("Example Neuron has less than 2 spikes, so it was omitted")
                    assert(1==2)
                end
            end
            cv(:, any(isnan(cv), 1)) = []; %omitnan
            if plot_name == "ex_neuron"
                plot(t(cv_window), cv(:, new_ex_neuron), 'Color', stim_color, 'Linewidth', 4)
            else
                plot(t(cv_window), mean(cv(:, 1:top_N), 2), 'Color', stim_color, 'Linewidth', 4)
            end
        end
        hold off
        ylabel("P1 CV")
        xlabel("Time (s)")
        legend(["Pulsatile", "Galvanic", "Control", "Anodic"])
    
    elseif plot_name == "p1_wins"
        figure(1);
        for sim_name = sim_names
            stim_cv = zeros(length(stim_amps), num_trials, num_group);
            load(sprintf("Simulation %s/ustim/r.mat", sim_name), 'ball_r')
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
                if contains(sim_name, "Discon")
                    stim_coherences = 0;
                end
                for trial = start_trial:end_trial
                    relative_trial = trial - start_trial + 1;
                    if ~contains(sim_name, "Discon") && (...
                            decisions(relative_trial, stim_coherences==c) ~= 1 || ...
                            decision_times(relative_trial, stim_coherences==c) > t_cut )
                        stim_cv(j, relative_trial, :) = NaN;
                        continue %skip trials where P1 doesn't win
                    end
                    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                        "recspikes")
                    isi = zeros(length(t), num_group);
                    for nn = 1:num_group
                        spikeidx = recspikes(int2str(nn));
                        if length(spikeidx) < 2
                            isi(:, nn) = NaN;
                            continue
                        end
                        neuron_isi_idx = diff(spikeidx);
                        pre_isi_idx = max(neuron_isi_idx(1), spikeidx(1));
                        isi(1:spikeidx(1)-1, nn) = t(pre_isi_idx);
                        for i = 1:length(neuron_isi_idx)-1
                            isi_idx = neuron_isi_idx(i);
                            spikei = spikeidx(i);
                            next_spikei = spikeidx(i+1);
                            isi(spikei:next_spikei-1, nn) = t(isi_idx);
                        end
                        post_isi_idx = max(neuron_isi_idx(end), length(t)-spikeidx(end));
                        isi(spikeidx(end):end, nn) = t(post_isi_idx);
                    end
                    cv = std(isi(cv_window, :), [], 1) ./ mean(isi(cv_window, :), 1);
                    stim_cv(j, relative_trial, :) = cv;
                end
            end
            pulse_cv = reshape(stim_cv(1, :, :), [num_trials, num_group]);
            galvanic_cv = reshape(stim_cv(2, :, :), [num_trials, num_group]);
            control_cv = reshape(stim_cv(3, :, :), [num_trials, num_group]);
            anodic_cv = reshape(stim_cv(4, :, :), [num_trials, num_group]);
            pulse_trialmean = mean(pulse_cv, 1, 'omitnan');
            galvanic_trialmean = mean(galvanic_cv, 1, 'omitnan');
            control_trialmean = mean(control_cv, 1, 'omitnan');
            anodic_trialmean = mean(anodic_cv, 1, 'omitnan');
            pulse_wins = all(~isnan(pulse_cv), 2);
            galvanic_wins = all(~isnan(galvanic_cv), 2);
            control_wins = all(~isnan(control_cv), 2);
            anodic_wins = all(~isnan(anodic_cv), 2);
            pulse_sem = std(pulse_cv, [], 1, 'omitnan') / sqrt(length(pulse_wins));
            galvanic_sem = std(galvanic_cv, [], 1, 'omitnan') / sqrt(length(galvanic_wins));
            control_sem =  std(control_cv, [], 1, 'omitnan') / sqrt(length(control_wins));
            anodic_sem = std(anodic_cv, [], 1, 'omitnan') / sqrt(length(anodic_wins));

            if contains(sim_name, "Discon")
                plot_shape = 'o';
            else
                plot_shape = '^';
            end
            figure(1);
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
            ylabel("Coefficient of Variation (unitless)")
        end
    end
end