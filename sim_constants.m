%%% Paul Adkisson
%%% 10/27/22
%%% Purpose: Define high-level simulation constants for each run

sim_name = "P1_Rec";
sim_path = sprintf("Simulation %s", sim_name);
sim_type = "p1_rec"; %"con" or "discon" or "ps_val" or "gs_val" or "p1_int" or "p1_rec"
start_trial = 29;
end_trial = 100;
pulse_coherences = [0] / 100;
control_coherences = [0] / 100;
galvanic_coherences = [0] / 100;
pulse_amps = [-10]*1e-6;
dc_amps = [-1.4, 0]*1e-6;
% dc_amps = [1.4]*1e-6;
%pulse_amps = [];
%dc_amps = []*1e-6;
bam_constants(sim_path, sim_type, start_trial, end_trial, pulse_coherences, ...
    control_coherences, galvanic_coherences, pulse_amps, dc_amps);
save(strcat(sim_path, "/sim_constants.mat"))
