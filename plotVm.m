sim_name = "Test";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))

load(strcat(sim_path, "/brain1/data/-10000.0nA_pulse/c=0.000/trial1.mat"), 'Vm')
figure;
plot(t(t>=1&t<1.01)*1e3, (Vm(t>=1&t<1.01, 1)-EL)*1e3)
xlabel("Time (ms)")
ylabel("Membrane Polarization Relative to V_{rest} (mV)")

figure;
plot(t, Vm(:, 1)*1e3)