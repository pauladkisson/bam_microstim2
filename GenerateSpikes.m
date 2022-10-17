%%% Paul Adkisson
%%% 9.2.21
%%% Purpose: Generate Poisson Spike Trains for various trials of a given
%%% frequency.
function GenerateSpikes(fr_bg, max_fr_task, coherences, f, N_E, N_I, t_task, ...
    t_taskoff, t, start_trial, end_trial, sim_path)
    dt = t(2) - t(1);
    N = N_E + N_I;
    num_group = floor(f*N_E);
    g1_idx = 1:num_group;
    g2_idx = num_group+1:num_group*2;
    time_idx = t>=t_task & t<t_taskoff;
    for c = coherences
        fr_task1 = max_fr_task/2*(1 + c);
        fr_task2 = max_fr_task/2*(1 - c);
        r0 = fr_bg;
        m = 0.5;
        f0 = 40;
        fr = ones(length(t), N) .* (r0 * (1 + m*cos(2*pi*f0*t'))); %Instaneous Firing Rate
        fr(time_idx, g1_idx) = fr(time_idx, g1_idx) + fr_task1;
        fr(time_idx, g2_idx) = fr(time_idx, g2_idx) + fr_task2;
        for trial = start_trial:end_trial
            rng(trial);
            spikes = rand(length(t), N) < (dt*fr);
            basepath = sprintf("/spikes/c=%0.3f/trial%0.0f", [c, trial]);
            spikepath = strcat(sim_path, basepath);
            mkdir(spikepath)
            save(strcat(spikepath, "/input.mat"), "spikes", "-v7.3")
        end
    end
end