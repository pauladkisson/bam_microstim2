function plot_psval(sim_path, plot_name, pulse_amp, ex_amp, ex_spont)
    load(strcat(sim_path, "/bam_constants.mat"))
    if plot_name == "amp"
        trial = find(fr_bgs==ex_spont); 
        load(strcat(sim_path, sprintf("/data/%0.2fuA_pulse/c=0.000/trial%0.0f.mat", [pulse_amp*1e6, trial])));
        lin_amps = true_amps(1)*1e6:-0.01:floor(true_amps(end)*1e6);
        spikes = zeros(length(t), N);
        for nn = 1:N
            for spike_idx = recspikes(int2str(nn))
                spikes(spike_idx, nn) = 1;
            end
        end
        neuron_frs = sum(spikes, 1) ./ t_span;
        cmap = flipud(parula(length(lin_amps))); %using the same colormap as thia
        figure;
        colormap(cmap)
        set(gca, 'fontsize', 18)
        hold on
        for i = 1:num_amps
            true_amp = true_amps(i);
            if i == 1
                start_idx = 1;
            else
                start_idx = floor(N/num_amps*(i-1)) + 1;
            end
            end_idx = floor(N/num_amps*i);
            plot(true_freqs, neuron_frs(start_idx:end_idx), 'Color', ...
                 cmap(abs(lin_amps-true_amp*1e6)<0.005, :), 'Linewidth', 2)
        end
        for i = 1:5
            plot(true_freqs, true_freqs/i, 'k--')
        end
        legend(compose("%0.1fuA", true_amps*1e6))
        xlim([0, 350])
        ylim([-20, 240])
        yticks(-20:20:240)
        hold off
        ylabel("Neuron Average Firing Rate (Hz)")
        xlabel("Pulse Stimulation Frequency (Hz)")
    elseif plot_name == "spont"
        delta_spont = 0.1;
        lin_sponts = 0:0.1:100;
        lin_cmap = turbo(length(lin_sponts));
        figure;
        colormap(lin_cmap)
        set(gca, 'fontsize', 18)
        hold on
        for trial = start_trial:end_trial
            load(strcat(sim_path, sprintf("/data/%0.2fuA_pulse/c=0.000/trial%0.0f.mat", [pulse_amp*1e6, trial])))
            spikes = zeros(length(t), N);
            for nn = 1:N
                for spike_idx = recspikes(int2str(nn))
                    spikes(spike_idx, nn) = 1;
                end
            end
            baseline_frs = sum(spikes(:, ps_stim_amps==0)) ./ t_span;
            neuron_frs = sum(spikes(:, abs(ps_stim_amps-ex_amp)<eps), 1) ./ t_span;
            delta_frs = neuron_frs - baseline_frs;
            plot(true_freqs, delta_frs, 'Color', lin_cmap(abs(lin_sponts-neuron_frs(1))<delta_spont/2, :))
        end
        for i = 1:5
            plot(true_freqs, true_freqs/i, 'k--')
        end
        hold off
        ylabel("Change in Neuron Average Firing Rate (Hz)")
        xlabel("Pulse Stimulation Frequency (Hz)")
        cb = colorbar('TickLabels', [lin_sponts(1), lin_sponts(end)], 'Direction', 'reverse');
        ylim([-100, 250])
        xlim([0, 350])
    end
end