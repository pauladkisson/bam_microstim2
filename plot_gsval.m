function plot_gsval(sim_path, galvanic_amp, plot_name)
    load(strcat(sim_path, "/bam_constants.mat"))
    delta_spont = 0.1;
    lin_sponts = 0:0.1:100;
    lin_cmap = turbo(length(lin_sponts));
    if plot_name == "internal"
        load(strcat(sim_path, sprintf("/ustim/%0.2fuA_galvanic.mat", galvanic_amp*1e6)))
        I_int = Vmir*gL(1);
    end
    figure;
    colormap(lin_cmap)
    set(gca, 'fontsize', 18)
    hold on
    for trial = start_trial:end_trial
        load(strcat(sim_path, sprintf("/data/%0.2fuA_galvanic/c=0.000/trial%0.0f.mat", [galvanic_amp*1e6, trial])))
        spikes = zeros(length(t), N);
        for nn = 1:N 
            for spike_idx = recspikes(int2str(nn))
                spikes(spike_idx, nn) = 1;
            end
        end
        baseline_frs = mean(sum(spikes(:, gs_stim_amps==0)) ./ t_span);
        neuron_frs = sum(spikes, 1) ./ t_span;
        delta_frs = neuron_frs - baseline_frs;
        avg_delta_frs = zeros(num_amps, 1);
        true_I_ints = zeros(num_amps, 1);
        for i = 1:num_amps
            start_idx = num_reps*(i-1) + 1;
            end_idx = num_reps*i;
            avg_delta_frs(i) = mean(delta_frs(start_idx:end_idx));
            if plot_name == "internal"
                true_I_ints(i) = I_int(start_idx);
            end
        end
        spont_color = lin_cmap(abs(lin_sponts-baseline_frs)<(delta_spont/2), :);
        if plot_name == "external"
            plot(true_amps*(-1e6), avg_delta_frs, 'Color', spont_color)
        else
            plot(true_I_ints*(-1e12), avg_delta_frs, 'Color', spont_color)
        end
    end
    hold off
    ylabel("Change in Neuron Average Firing Rate (Hz)")
    if plot_name == "external"
        xlabel("Cathodic Amplitude (uA)")
    else
        xlabel("Internal Stimulation Amplitude (pA)")
    end
    cb = colorbar('TickLabels', [lin_sponts(1), lin_sponts(end)], 'Direction', 'reverse');
    ylim([-100, 250])
end