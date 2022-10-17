%%% Paul Adkisson
%%% 6.24.22
%%% Copies files from Rita to external hard drive
%%% Note: Pulse amps and Galvanic amps are in nA, unlike in the main
%%% scripts.
function rita2ext(root_source, root_dest, sim_name, pulse_amps, galvanic_amps, brains)
    sim_source = fullfile(root_source, sprintf("Simulation %s", sim_name));
    sim_dest = fullfile(root_dest, sprintf("Simulation %s", sim_name));
    for brain = brains
        brain_source = fullfile(sim_source, sprintf("brain%0.0f", brain));
        brain_dest = fullfile(sim_dest, sprintf("brain%0.0f", brain));
        for p_amp = pulse_amps
           p_source = fullfile(brain_source, "data", sprintf("%0.1fnA_pulse", p_amp));
           p_dest = fullfile(brain_dest, "data", sprintf("%0.1fnA_pulse", p_amp));
           [status, msg] = copyfile(p_source, p_dest);
           assert(status, msg)
        end
        for g_amp = galvanic_amps
           g_source = fullfile(brain_source, "data", sprintf("%0.1fnA_galvanic", g_amp));
           g_dest = fullfile(brain_dest, "data", sprintf("%0.1fnA_galvanic", g_amp));
           [status, msg] = copyfile(g_source, g_dest);
           assert(status, msg)
        end
    end
end