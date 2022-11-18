%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Synchrony
function plot_sync(sim_names, pulse_amps, stim_amps, t, t_cut, num_group, num_affected, ...
                        N_start, N_end, win_start, win_stop, c_win, ex_c, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        anodic_coherences, start_trial, end_trial, num_trials, symmetric)
    for sim_name = sim_names
        disp(sim_name)
        
        stim_sync = zeros(length(stim_amps), num_trials, num_group, num_group);
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
                fprintf("Trial %0.0f \n", trial)
                if ~contains(sim_name, "Discon") && (...
                            decisions(relative_trial, stim_coherences==c) ~= 1 || ...
                            decision_times(relative_trial, stim_coherences==c) > t_cut )
                    stim_sync(j, trial, :, :) = NaN;
                    continue %skip trials where P1 doesn't win or decision takes too long
                end
                load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                    "recspikes")
                if symmetric
                    pairwise_sync = get_sym_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop);
                else
                    pairwise_sync = get_pairwise_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop);
                end
                stim_sync(j, trial, :, :) = pairwise_sync;
            end
        end
        pulse_sync = reshape(mean(stim_sync(1, :, :, :), 2, 'omitnan'), [num_group, num_group]);
        galvanic_sync = reshape(mean(stim_sync(2, :, :, :), 2, 'omitnan'), [num_group, num_group]);
        control_sync = reshape(mean(stim_sync(3, :, :, :), 2, 'omitnan'), [num_group, num_group]);
        anodic_sync = reshape(mean(stim_sync(4, :, :, :), 2, 'omitnan'), [num_group, num_group]);

        nan_color = uint8([0, 0, 128]);
        ticks = log10([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]);
        tick_labels = ["1", "10", "", "", "", "", "", "", "", "", "100"];
        figure;
        h = heatmap(log10(pulse_sync(1:num_affected, 1:num_affected)*100), ...
                    'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
        colormap(hot)
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Pulsatile")
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = ticks;
        cb.TickLabels = tick_labels;

        figure;
        h = heatmap(log10(galvanic_sync(1:num_affected, 1:num_affected)*100), ...
                    'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
        colormap(hot)
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Galvanic")
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = ticks;
        cb.TickLabels = tick_labels;

        figure;
        h = heatmap(log10(control_sync(1:num_affected, 1:num_affected)*100), ...
                    'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
        colormap(hot)
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Control")
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = ticks;
        cb.TickLabels = tick_labels;
        
        figure;
        h = heatmap(log10(anodic_sync(1:num_affected, 1:num_affected)*100), ...
                    'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
        colormap(hot)
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Anodic")
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = ticks;
        cb.TickLabels = tick_labels;
    end
end