function [frs, fr_var] = recspikes2popfrs(recspikes, t, N, dt, p, f, N_E)
    spikes = zeros(length(t), N);
    for nn = 1:N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    [frs, fr_var] = spikes2popfrs(spikes, dt, p, f, N_E);
end