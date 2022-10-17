%%% Paul Adkisson
%%% 6.24.22
%%% Compresses data by removing pop_frs
sim_name = "Test";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))


start_trial = 1;
end_trial = 4;
brains = 1;
%pulse_coherences = [-100, -78.8, -75.6, -72.4, -69.2, -66, -51.2, -25.6, 0, 25.6] / 100;
%control_coherences = [-100, -51.2, -25.6, -12.8, -6.4, -3.2, 0, 3.2, 6.4, 12.8, 25.6] / 100;
%galvanic_coherences = [-100, -51.2 -42.6, -39.4, -36.2, -33, -29.8, -25.6, 0, 25.6] / 100;
%galvanic_coherences = [-100, -75.6, -51.2, -25.6, 0, 25.6] / 100;
pulse_coherences = [0];
control_coherences = [0];
galvanic_coherences = [0];
pulse_amps = [];
dc_amps = [-28]*1e-9;
stim_amps = [pulse_amps, dc_amps];
%}

for brain = brains
    fprintf("Brain %0.0f \n", brain)
    for k = 1:length(stim_amps)
        stim_amp = stim_amps(k);
        pulse = k <= length(pulse_amps);
        if pulse
            fprintf("Pulse Stimulation Amplitude: %0.1fnA \n", stim_amp*1e9)
            output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                [sim_name, brain, stim_amp*1e9]);
            coherences = pulse_coherences;
        else
            fprintf("Galvanic Stimulation Amplitude: %0.1fnA \n", stim_amp*1e9)
            output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                [sim_name, brain, stim_amp*1e9]);
            if stim_amp == 0
                coherences = control_coherences;
                if brain ~= 1
                    continue
                end
            else
                coherences = galvanic_coherences;
            end
        end
        for i = 1:length(coherences)
            c = coherences(i);
            fprintf("Coherence: %0.1f%% \n", c*100)
            output_coherentpath = strcat(output_stimpath, sprintf("/c=%0.3f", c));
            for trial = start_trial:end_trial
                output_trialpath = strcat(output_coherentpath, sprintf("/trial%0.0f.mat", trial));
                load(output_trialpath, "recspikes")
                save(output_trialpath, "recspikes")
            end
        end
    end
end