function plot_gsval(sim_path, galvanic_amp, plot_name)
    load(strcat(sim_path, "/bam_constants.mat"))
    cmap = turbo(end_trial);
    lin_cmap = turbo(1000);
    if plot_name == "internal"
        load(strcat(sim_path, sprintf("/ustim/%0.2fuA_galvanic.mat", galvanic_amp*1e6)))
        I_int = Vmir*gL(1);
    end
    figure;
    colormap(lin_cmap)
    set(gca, 'fontsize', 18)
    hold on
    spont_fr_lims = [0, 0];
    for trial = start_trial:end_trial
        load(strcat(sim_path, sprintf("/data/%0.2fuA_galvanic/c=0.000/trial%0.0f.mat", [galvanic_amp*1e6, trial])))
        spikes = zeros(length(t), N);
        for nn = 1:N 
            for spike_idx = recspikes(int2str(nn))
                spikes(spike_idx, nn) = 1;
            end
        end
        baseline_frs = sum(spikes(:, true_amps==0)) ./ t_span;
        neuron_frs = sum(spikes, 1) ./ t_span;
        delta_frs = neuron_frs - baseline_frs;
        avg_delta_frs = zeros(num_amps, 1);
        true_I_ints = zeros(num_amps, 1);
        for i = 1:num_amps
            start_idx = num_reps*(i-1) + 1;
            end_idx = num_reps*i;
            avg_delta_frs(i) = mean(delta_frs(start_idx:end_idx));
            true_I_ints(i) = I_int(start_idx);
        end

        if plot_name == "external"
            plot(abs(true_amps*1e6), avg_delta_frs, 'Color', cmap(trial, :))
        else
            plot(true_I_ints*1e12, avg_delta_frs, 'Color', cmap(trial, :))
        end
        if trial == start_trial
            spont_fr_lims(1) = neuron_frs(1);
        elseif trial == end_trial
            spont_fr_lims(2) = neuron_frs(1);
        end
    end
    hold off
    ylabel("Change in Neuron Average Firing Rate (Hz)")
    if plot_name == "external"
        xlabel("Cathodic Amplitude (uA)")
    else
        xlabel("Internal Stimulation Amplitude (pA)")
    end
    cb = colorbar('TickLabels', [spont_fr_lims(1), spont_fr_lims(end)], 'Direction', 'reverse');
end