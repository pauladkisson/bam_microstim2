%%% Paul Adkisson
%%% 9.2.21
%%% Purpose: Generate Poisson Spike Trains for various trials of a vector 
%%% of frequencies.
function GeneratePopSpikes(sim_path, fr_bgs, coherences, N, t, ...
    start_trial, end_trial)
    dt = t(2) - t(1);
    for c = coherences
        coherentpath = strcat(sim_path, sprintf("/spikes/c=%0.3f", c));
        mkdir(coherentpath);
        for trial = start_trial:end_trial
            fr_bg = fr_bgs(trial);
            spikes = rand(length(t), N) < (dt*fr_bg);
            spikepath = strcat(coherentpath, sprintf("/trial%0.0f.mat", trial));
            save(spikepath, 'spikes', '-v7.3');
        end
    end
end