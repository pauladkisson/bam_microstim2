%%% Paul Adkisson
%%% 9.2.21
%%% Purpose: Generate Poisson Spike Trains for various trials of a given
%%% frequency.
function GenerateSpikes(fr_bg, m, f0, max_fr_task, coherences, f, p, N_E, N_I, ...
        t_task, t_taskoff, t, start_trial, end_trial, sim_path)
    dt = t(2) - t(1);
    N = N_E + N_I;
    num_group = floor(f*N_E);
    g1_idx = 1:num_group;
    g2_idx = num_group+1:num_group*p;
    time_idx = t>=t_task & t<t_taskoff;
    for c = coherences
        fr_task1 = max_fr_task/2*(1 + c);
        fr_task2 = max_fr_task/2*(1 - c);
        r0 = fr_bg;
        fr = ones(length(t), N) .* (r0 * (1 + m*cos(2*pi*f0*t'))); %Instaneous Firing Rate
        fr(time_idx, g1_idx) = fr(time_idx, g1_idx) + fr_task1;
        fr(time_idx, g2_idx) = fr(time_idx, g2_idx) + fr_task2;
        coherentpath = strcat(sim_path, sprintf("/spikes/c=%0.3f", c));
        mkdir(coherentpath);
        for trial = start_trial:end_trial
            rng(trial);
            spikes = rand(length(t), N) < (dt*fr);
            spikepath = strcat(coherentpath, sprintf("/trial%0.0f.mat", trial));
            save(spikepath, 'spikes', '-v7.3');
        end
    end
end