%Identical to get pairwise sync except symmetric where only the fraction of
%spikes from the neuron with more spikes is reported
function sym_sync = get_sym_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop)
    N = N_end - N_start + 1;
    sym_sync = zeros(N, N);
    sym_sync(logical(eye(N, N))) = NaN;
    for i = N_start:N_end
        idx_i = i - N_start + 1;
        n1_spiketimes = t(recspikes(int2str(i)));
        n1_spiketimes = n1_spiketimes(n1_spiketimes>=win_start & n1_spiketimes<win_stop);
        for j = i+1:N_end
            idx_j = j - N_start + 1;
            n2_spiketimes = t(recspikes(int2str(j)));
            n2_spiketimes = n2_spiketimes(n2_spiketimes>=win_start & n2_spiketimes<win_stop);
            if length(n1_spiketimes) >= length(n2_spiketimes) 
                [sync, ~] = get_sync(n1_spiketimes, n2_spiketimes, c_win);
            else
                [sync, ~] = get_sync(n2_spiketimes, n1_spiketimes, c_win);
            end
            sym_sync(idx_i, idx_j) = sync;
            sym_sync(idx_j, idx_i) = sync;
        end
    end
    sym_sync(isnan(sym_sync')) = NaN;
end