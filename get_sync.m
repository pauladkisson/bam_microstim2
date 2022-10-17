function [sync, coincident_spikes] = get_sync(n1_spiketimes, n2_spiketimes, c_win)
    coincident_spikes = zeros(size(n1_spiketimes));
    for i = 1:length(n1_spiketimes)
        n1_spiketime = n1_spiketimes(i);
        t_diff = min(abs(n1_spiketime - n2_spiketimes));
        if t_diff < c_win
            coincident_spikes(i) = 1;
        end
    end
    sync = sum(coincident_spikes) / length(n1_spiketimes);
end