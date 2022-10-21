function [frs, fr_var] = spikes2popfrs(spikes, dt, p, f, N_E)
    win_size = 5e-3;
    avg_win_size = 50e-3;
    num_group = floor(f*N_E);
    w = ones(floor(win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, spikes) ./ dt;
    w = ones(floor(avg_win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, neuron_frs);
    frs = zeros(size(spikes, 1), p+2);
    fr_var = zeros(size(spikes, 1), p+2);
    for i = 1:p
        group_idx = ((i-1)*num_group+1):i*num_group;
        frs(:, i) = mean(neuron_frs(:, group_idx), 2);
        fr_var(:, i) = std(neuron_frs(:, group_idx), [], 2);
    end
    frs(:, end-1) = mean(neuron_frs(:, f*p*N_E+1:N_E), 2);
    fr_var(:, end-1) = std(neuron_frs(:, f*p*N_E+1:N_E), [], 2);
    frs(:, end) = mean(neuron_frs(:, N_E+1:end), 2);
    fr_var(:, end) = std(neuron_frs(:, N_E+1:end), [], 2);
end