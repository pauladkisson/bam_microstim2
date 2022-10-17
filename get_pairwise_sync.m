function pairwise_sync = get_pairwise_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop)
    N = N_end - N_start + 1;
    pairwise_sync = zeros(N, N);
    pairwise_sync(logical(eye(N, N))) = NaN;
    for i = N_start:N_end
        idx_i = i - N_start + 1;
        n1_spiketimes = t(recspikes(int2str(i)));
        n1_spiketimes = n1_spiketimes(n1_spiketimes>=win_start & n1_spiketimes<win_stop);
        for j = N_start:N_end
            idx_j = j - N_start + 1;
            if i == j
                continue
            end
            n2_spiketimes = t(recspikes(int2str(j)));
            n2_spiketimes = n2_spiketimes(n2_spiketimes>=win_start & n2_spiketimes<win_stop);
            [sync, ~] = get_sync(n1_spiketimes, n2_spiketimes, c_win);
            pairwise_sync(idx_i, idx_j) = sync;
        end
    end
    pairwise_sync(isnan(pairwise_sync')) = NaN;
end