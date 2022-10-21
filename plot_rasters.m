%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Example Rasters

function plot_rasters(sim_name, pulse_amps, stim_amps, ex_neurons, t, ...
                      t_task, t_taskoff, stim_freq, default_colors, top_N, ...
                      ex_stim_j, ex_c, ex_trial, tlim, plot_name)
    load(sprintf("Simulation %s/ustim/r.mat", sim_name), "ball_r")
    if plot_name == "subplot"
        figure;
    end
    for j = 1:length(stim_amps)
        if plot_name == "single_stim" && ex_stim_j ~= j || plot_name == "grouped_stim" 
            continue
        end
        stim_amp = stim_amps(j);
        pulse = j<=length(pulse_amps);
        if pulse
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                [sim_name, stim_amp*1e6]);
        else
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e6]);
        end
        load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), ...
            "recspikes")
        spikes = zeros(length(t), top_N);
        for nn = 1:top_N
            for spike_idx = recspikes(int2str(nn))
                spikes(spike_idx, nn) = 1;
            end
        end

        [g1_time_idx, g1_neuron_idx, g1_idx] = get_spike_idx(spikes);
        
        if plot_name == "subplot"
            subplot(3, 1, j)
            hold on
        else
            figure;
            hold on
        end
        scatter(t(g1_time_idx), ball_r(g1_idx(g1_neuron_idx))*1e6, 1, "Marker", "|", ...
            "MarkerFaceColor", default_colors(1, :), "MarkerEdgeColor", default_colors(1, :))
        xlim(tlim)
        xlabel("Time (s)")
        ylabel("Distance from Electrode (um)")
        if pulse
            title("Pulsatile Stimulation")
            pulsetimes = t_task:1/stim_freq:(t_taskoff-1/stim_freq);
            scatter(pulsetimes, -5, "k", "Marker", "|")
            xlim(tlim)
        else
            if abs(stim_amp) > 0
                title("Galvanic Stimulation")
            else
                title("Control")
            end
        end
    end
    set(findall(gcf,'-property','FontSize'),'FontSize', 18)
    hold off
    
    if plot_name == "grouped_stim"
        top_N = length(ex_neurons);
        neuron_pos = 1:length(ex_neurons);
        raster_size = 100;
        pulse_color = [0, 0, 1];
        figure;
        set(gca, 'fontsize', 20)
        hold on
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            pulse = j<=length(pulse_amps);
            if pulse
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                    [sim_name, stim_amp*1e6]);
                spike_color = default_colors(7, :);
            else
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                    [sim_name, stim_amp*1e6]);
                if stim_amp == 0
                    spike_color = [0, 0, 0];
                else
                    spike_color = default_colors(5, :);
                end
            end
            load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), ...
                "recspikes")
            spikes = zeros(length(t), top_N);
            for nn = 1:top_N
                ex_neuron = ex_neurons(nn);
                for spike_idx = recspikes(int2str(ex_neuron))
                    spikes(spike_idx, nn) = 1;
                end
            end
            [g1_time_idx, g1_neuron_idx, ~] = get_spike_idx(spikes);
            pos = neuron_pos(g1_neuron_idx);
            scatter(t(g1_time_idx), (pos-1)*5+j, raster_size, "Marker", "|", ...
                "MarkerFaceColor", spike_color, "MarkerEdgeColor", spike_color)
        end

        pulsetimes = t_task:1/stim_freq:(t_taskoff-1/stim_freq);
        pulse_pos = ones(length(pulsetimes), length(neuron_pos));
        for nn = 1:length(neuron_pos)
            pulse_pos(:, nn) = (neuron_pos(nn)-1)*5;
        end
        scatter(pulsetimes, pulse_pos, raster_size, "Marker", "|", ...
            "MarkerFaceColor", pulse_color, "MarkerEdgeColor", pulse_color)
        hold off
        xlabel("Time (s)")
        xlim(tlim)
        yticks([])
        yticklabels([])
    end
end

function [time_idx, neuron_idx, g_idx] = get_spike_idx(g_spikes)
    g_idx = 1:size(g_spikes, 2);
    g_spikes = g_spikes(:, g_idx);
    [time_idx, neuron_idx] = find(g_spikes);
end