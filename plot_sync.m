%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Synchrony
function plot_sync(sim_names, pulse_amps, stim_amps, t, num_group, N_start, ...
                        N_end, win_start, win_stop, c_win, c, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        start_trial, end_trial, num_trials, symmetric)
    for sim_name = sim_names
        disp(sim_name)
        
        stim_sync = zeros(length(stim_amps), num_trials, num_group, num_group);
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            pulse = j<=length(pulse_amps);
            if pulse
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                    [sim_name stim_amp*1e6]);
                %stim_coherences = pulse_coherences;
                disp("Pulsatile")
            else
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                    [sim_name, stim_amp*1e6]);
                if stim_amp == 0
                    %stim_coherences = control_coherences;
                    disp("Control")
                else
                    %stim_coherences = galvanic_coherences;
                    disp("Galvanic")
                end
            end
            %load(strcat(output_stimpath, "/decisions.mat"), "decisions")
            load(sprintf("Simulation %s/data/0.00uA_galvanic/decisions.mat", ...
                    sim_name), 'decisions'); %using only control decisions
            stim_coherences = control_coherences;
            for trial = start_trial:end_trial
                fprintf("Trial %0.0f \n", trial)
                if sim_name==sim_names(1) && decisions(trial, stim_coherences==c) ~= 1
                    disp("P1 lost")
                    stim_sync(j, trial, :, :) = NaN;
                    continue %skip trials where P1 doesn't win
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

        nan_color = uint8([0, 0, 128]);
        ticks = log10([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]);
        tick_labels = ["1", "10", "", "", "", "", "", "", "", "", "100"];
        figure;
        h = heatmap(log10(pulse_sync*100), 'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
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
        h = heatmap(log10(galvanic_sync*100), 'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
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
        h = heatmap(log10(control_sync*100), 'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
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
    end
end